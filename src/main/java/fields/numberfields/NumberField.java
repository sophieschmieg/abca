package fields.numberfields;

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.SortedMap;
import java.util.SortedSet;
import java.util.TreeMap;
import java.util.TreeSet;

import fields.finitefields.FiniteField;
import fields.finitefields.FiniteField.FFE;
import fields.finitefields.PrimeField;
import fields.finitefields.PrimeField.PFE;
import fields.floatingpoint.Complex;
import fields.floatingpoint.Complex.ComplexNumber;
import fields.floatingpoint.FiniteRealVectorSpace;
import fields.floatingpoint.Reals;
import fields.floatingpoint.Reals.Real;
import fields.helper.AbstractElement;
import fields.helper.AbstractFieldExtension;
import fields.helper.FieldEmbedding;
import fields.helper.GaloisGroup;
import fields.helper.GenericAlgebraicRingExtension.GenericAlgebraicExtensionElement;
import fields.integers.Integers;
import fields.integers.Integers.IntE;
import fields.integers.Rationals;
import fields.integers.Rationals.Fraction;
import fields.integers.ValueFractions;
import fields.interfaces.AlgebraicExtensionElement;
import fields.interfaces.DiscreteValuationField.OtherVersion;
import fields.interfaces.UnivariatePolynomialRing.ExtendedResultantResult;
import fields.interfaces.DiscreteValuationRing;
import fields.interfaces.GlobalFieldExtension;
import fields.interfaces.Ideal;
import fields.interfaces.MathMap;
import fields.interfaces.Polynomial;
import fields.interfaces.PolynomialRing;
import fields.interfaces.UnivariatePolynomial;
import fields.interfaces.UnivariatePolynomialRing;
import fields.local.PAdicField;
import fields.local.Value;
import fields.numberfields.CompletedNumberField.Ext;
import fields.numberfields.IdealClassGroup.IdealClass;
import fields.numberfields.NumberField.NFE;
import fields.numberfields.NumberFieldIntegers.NumberFieldIdeal;
import fields.polynomials.AbstractPolynomialRing;
import fields.polynomials.Monomial;
import fields.vectors.Matrix;
import fields.vectors.RealLattice;
import fields.vectors.Vector;
import util.FunctionMathMap;
import util.MiscAlgorithms;
import util.Pair;
import util.SingletonSortedMap;
import varieties.curves.elliptic.EllipticCurve;

public class NumberField extends AbstractFieldExtension<Fraction, NFE, NumberField> implements
		GlobalFieldExtension<Fraction, IntE, PFE, NFE, NFE, PFE, FFE, FiniteField, LocalizedNumberField, NumberField> {
	private static final int PRECISION_FOR_EMBEDDINGS = 128;

	private Rationals q;
	private Reals r;
	private Complex c;
	private boolean embeddingsComputed;
	private List<EmbeddedNumberField<Real, Reals>> realEmbeddings;
	private List<EmbeddedNumberField<ComplexNumber, Complex>> complexEmbeddings;
	private Map<IntE, Map<Integer, List<EmbeddedNumberField<Ext, CompletedNumberField>>>> padicEmbeddings;
	private FiniteRealVectorSpace logRepresentationSpace;
	private FiniteRealVectorSpace minkowskiEmbeddingSpace;
	private NumberFieldIntegers maximalOrder;
	private IdealGroup idealGroup;
	private IdealClassGroup idealClassGroup;
	private Map<Polynomial<NFE>, List<Polynomial<NFE>>> polynomialFactorizationCache;
	private NumberFieldIdeal idealForFactorization;
	private boolean hasIntegerMinimalPolynomial;
	private NumberField withIntegerMinimalPolynomial;
	private IntE multiplier;
	private static Map<Polynomial<Fraction>, NumberField> numberFields = new TreeMap<>();
	private GaloisGroup<Fraction, NFE, NumberField> galoisGroup;

	public static class NFE extends AbstractElement<NFE> implements AlgebraicExtensionElement<Fraction, NFE> {
		private UnivariatePolynomial<Fraction> e;

		NFE(NumberField field, UnivariatePolynomial<Fraction> e) {
			this.e = e;
		}

		@Override
		public int compareTo(NFE o) {
			return e.compareTo(((NFE) o).e);
		}

		@Override
		public String toString() {
			return this.e.toString(/* "α", */true);
		}

		@Override
		public UnivariatePolynomial<Fraction> asPolynomial() {
			return e;
		}

	}

	public static NumberField getNumberField(UnivariatePolynomial<Fraction> minimalPolynomial) {
		if (!minimalPolynomial.leadingCoefficient().equals(Rationals.q().one())) {
			throw new ArithmeticException("Non normal minimal polynomial!");
		}
		if (!numberFields.containsKey(minimalPolynomial)) {
			numberFields.put(minimalPolynomial, new NumberField(minimalPolynomial));
		}
		return numberFields.get(minimalPolynomial);
	}

	public static NumberField getNumberField() {
		UnivariatePolynomial<Fraction> minimalPolynomial = Rationals.q().getUnivariatePolynomialRing().getVar();
		if (!numberFields.containsKey(minimalPolynomial)) {
			numberFields.put(minimalPolynomial, new NumberField());
		}
		return numberFields.get(minimalPolynomial);
	}

	public static NumberField getNumberFieldFromIntegerPolynomial(UnivariatePolynomial<IntE> minimalPolynomial) {
		Rationals q = Rationals.q();
		UnivariatePolynomialRing<Fraction> polynomialRing = q.getUnivariatePolynomialRing();
		return getNumberField(
				polynomialRing.normalize(polynomialRing.getEmbedding(minimalPolynomial, q.getEmbeddingMap())));
	}

	public static NumberField getRootsOfUnityField(int order) {
		Integers z = Integers.z();
		FactorizationResult<IntE, IntE> factors = z.uniqueFactorization(z.getInteger(order));
		NumberField result = getNumberField();
		for (IntE prime : factors.primeFactors()) {
			UnivariatePolynomialRing<NFE> polynomialRing = result.getUnivariatePolynomialRing();
			UnivariatePolynomial<NFE> primeUnityPolynomial = polynomialRing.toUnivariate(polynomialRing.divideChecked(
					polynomialRing.subtract(polynomialRing.getVarPower(prime.intValueExact()), polynomialRing.one()),
					polynomialRing.subtract(polynomialRing.getVar(), polynomialRing.one())));
			FieldEmbedding<Fraction, NFE, NumberField> primeUnity = result.getEmbeddedExtension(primeUnityPolynomial);
			result = primeUnity.getField();
			polynomialRing = result.getUnivariatePolynomialRing();
			int remainingPower = z.power(prime, factors.multiplicity(prime) - 1).intValueExact();
			UnivariatePolynomial<NFE> primePowerUnity = polynomialRing
					.toUnivariate(polynomialRing.subtract(polynomialRing.getVarPower(remainingPower),
							polynomialRing.getEmbedding(primeUnity.getGenerator())));
			result = result.getEmbeddedExtension(primePowerUnity).getField();
		}
		return result;
	}

	public static NumberField getRealRootsOfUnityField(int order) {
		NumberField rootsOfUnity = getRootsOfUnityField(order);
		NFE primitiveRoot = rootsOfUnity.maximalOrder().primitiveRootsOfUnity().get(Integers.z().getInteger(order));
		NFE real = rootsOfUnity.add(primitiveRoot, rootsOfUnity.inverse(primitiveRoot));
		return getNumberField(rootsOfUnity.minimalPolynomial(real));
	}

	private NumberField() {
		super(Rationals.q());
		init();
	}

	private static String variableName(UnivariatePolynomial<Fraction> minimalPolynomial) {
		if (minimalPolynomial.degree() != 2) {
			return "x";
		}
		Rationals q = Rationals.q();
		if (minimalPolynomial.univariateCoefficient(0).equals(q.one())
				&& minimalPolynomial.univariateCoefficient(1).equals(q.zero())
				&& minimalPolynomial.univariateCoefficient(2).equals(q.one())) {
			return "i";
		}
		return "x";
	}

	private NumberField(UnivariatePolynomial<Fraction> minimalPolynomial) {
		super(minimalPolynomial, Rationals.q(), variableName(minimalPolynomial));
		init();
	}

	private void init() {
		this.q = Rationals.q();
		this.r = Reals.r(PRECISION_FOR_EMBEDDINGS);
		this.c = Complex.c(PRECISION_FOR_EMBEDDINGS);
		this.realEmbeddings = new ArrayList<>();
		this.complexEmbeddings = new ArrayList<>();
		this.embeddingsComputed = false;
		this.polynomialFactorizationCache = new TreeMap<>();
		Integers z = Integers.z();
		hasIntegerMinimalPolynomial = true;
		IntE lcm = z.one();
		for (int i = 0; i < degree(); i++) {
			Fraction coeff = minimalPolynomial().univariateCoefficient(i);
			if (!q.isInteger(coeff)) {
				hasIntegerMinimalPolynomial = false;
				IntE denominator = coeff.getDenominator();
				FactorizationResult<IntE, IntE> denominatorFactors = z.uniqueFactorization(denominator);
				for (IntE prime : denominatorFactors.primeFactors()) {
					lcm = z.lcm(
							z.power(prime,
									MiscAlgorithms.DivRoundUp(denominatorFactors.multiplicity(prime), degree() - i)),
							lcm);
				}
			}
		}
		if (!hasIntegerMinimalPolynomial) {
			UnivariatePolynomial<Fraction> integerMinimalPolynomial = q.getUnivariatePolynomialRing()
					.substitute(minimalPolynomial(), Collections.singletonList(
							q.getUnivariatePolynomialRing().getEmbedding(q.inverse(q.getInteger(lcm)), 1)));
			integerMinimalPolynomial = q.getUnivariatePolynomialRing().normalize(integerMinimalPolynomial);
			withIntegerMinimalPolynomial = getNumberField(integerMinimalPolynomial);
			multiplier = lcm;
		}
//		if (!q.isIrreducible(minimalPolynomial())) {
//			throw new ArithmeticException("Not irreducible");
//		}
	}

	public boolean hasIntegerMinimalPolynomial() {
		return hasIntegerMinimalPolynomial;
	}

	public NumberField withIntegerMinimalPolynomial() {
		return hasIntegerMinimalPolynomial ? this : withIntegerMinimalPolynomial;
	}

	public NFE toIntegerMinimalPolynomial(NFE t) {
		return hasIntegerMinimalPolynomial ? t : divide(t, getInteger(multiplier));
	}

	public NFE fromIntegerMinimalPolynomial(NFE t) {
		return hasIntegerMinimalPolynomial ? t : multiply(multiplier, t);
	}

	public MathMap<NFE, NFE> toIntegerMinimalPolynomialMap() {
		return new FunctionMathMap<>((NFE t) -> toIntegerMinimalPolynomial(t));
	}

	public MathMap<NFE, NFE> fromIntegerMinimalPolynomialMap() {
		return new FunctionMathMap<>((NFE t) -> fromIntegerMinimalPolynomial(t));
	}

	@Override
	public NumberFieldIntegers ringOfIntegers() {
		return maximalOrder();
	}

	@Override
	public Rationals getBaseField() {
		return q;
	}

	@Override
	public NFE getInteger(NFE t) {
		return t;
	}

	public NumberFieldIntegers maximalOrder() {
		if (maximalOrder == null) {
			if (!hasIntegerMinimalPolynomial()) {
				throw new ArithmeticException("Does not have integer minimal polynomial!");
			}
			maximalOrder = new NumberFieldIntegers(this);
		}
		return maximalOrder;
	}

	public NumberFieldOrder getOrder(NFE generator) {
		return getOrder(Collections.singletonList(generator));
	}

	public NumberFieldOrder getOrder(List<NFE> generators) {
		return new NumberFieldOrder(this, generators);
	}

	@Override
	public NumberField makeExtension(UnivariatePolynomial<Fraction> minimalPolynomial) {
		return getNumberField(minimalPolynomial);
	}

	@Override
	public ExtensionOfGlobalField<NFE, NFE, FFE, Fraction, IntE, PFE, NFE, NFE, PFE, FFE, FiniteField, LocalizedNumberField, NumberField> getGlobalFieldExtension(
			UnivariatePolynomial<NFE> minimalPolynomial) {
		Extension<NFE, Fraction, NFE, NumberField> extension = getExtension(minimalPolynomial);
		return new ExtensionOfGlobalField<>(this, extension.extension(), extension.embeddingMap(),
				extension.asVectorMap());
	}

	@Override
	public Exactness exactness() {
		return Exactness.EXACT;
	}

	public NFE getEmbedding(IntE t) {
		return getEmbedding(q.getEmbedding(t));
	}

	public NFE getEmbedding(BigInteger t) {
		return getEmbedding(Integers.z().getInteger(t));
	}

	@Override
	public NFE getInteger(BigInteger t) {
		return getEmbedding(t);
	}

	public boolean isInteger(NFE t) {
		return maximalOrder().isElement(t);
	}

	private void computeEmbeddings() {
		if (embeddingsComputed) {
			return;
		}
		Polynomial<ComplexNumber> minPoly = MiscAlgorithms.mapPolynomial(minimalPolynomial(), new MathMap<>() {
			@Override
			public ComplexNumber evaluate(Fraction t) {
				return c.getEmbedding(r.getEmbedding(t));
			}
		}, c.getUnivariatePolynomialRing());
		Map<ComplexNumber, Integer> roots = c.roots(minPoly);
		for (ComplexNumber alpha : roots.keySet()) {
			if (alpha.complexPart().equals(r.zero()) || alpha.complexPart().exponent() < -r.precision() + 1) {
				MathMap<Real, NFE> round = degree() == 1 ? new MathMap<>() {
					@Override
					public NFE evaluate(Real t) {
						return getInteger(t.round());
					}
				} : null;
				realEmbeddings.add(new EmbeddedNumberField<>(this, r, alpha.realPart(), round));
			} else if (alpha.complexPart().compareTo(r.zero()) > 0) {
				MathMap<ComplexNumber, NFE> round = degree() == 2 ? new MathMap<>() {
					@Override
					public NFE evaluate(ComplexNumber t) {
						IntE alphaPart = r.divide(t.complexPart(), alpha.complexPart()).round();
						ComplexNumber residue = c.subtract(t, c.multiply(alphaPart, alpha));
						IntE realPart = residue.realPart().round();
						return add(getInteger(realPart), multiply(alphaPart, alpha()));
					}
				} : null;
				complexEmbeddings.add(new EmbeddedNumberField<>(this, c, alpha, round));
			}
		}
		embeddingsComputed = true;
		if (realEmbeddings.size() + 2 * complexEmbeddings.size() != degree()) {
			throw new ArithmeticException("Embeddings not computed correctly!");
		}
	}

	public List<EmbeddedNumberField<Real, Reals>> realEmbeddings() {
		computeEmbeddings();
		return Collections.unmodifiableList(this.realEmbeddings);
	}

	public List<EmbeddedNumberField<ComplexNumber, Complex>> complexEmbeddings() {
		computeEmbeddings();
		return Collections.unmodifiableList(this.complexEmbeddings);
	}

	public List<EmbeddedNumberField<ComplexNumber, Complex>> complexEmbeddings(Complex c) {
		List<EmbeddedNumberField<ComplexNumber, Complex>> complexEmbeddings = new ArrayList<>();
		Polynomial<ComplexNumber> minPoly = MiscAlgorithms.mapPolynomial(minimalPolynomial(), new MathMap<>() {
			@Override
			public ComplexNumber evaluate(Fraction t) {
				return c.getEmbedding(r.getEmbedding(t));
			}
		}, c.getUnivariatePolynomialRing());
		Map<ComplexNumber, Integer> roots = c.roots(minPoly);
		for (ComplexNumber alpha : roots.keySet()) {
			if (alpha.complexPart().equals(r.zero()) || alpha.complexPart().exponent() < -r.precision() + 5) {
				continue;
			} else if (alpha.complexPart().compareTo(r.zero()) > 0) {
				MathMap<ComplexNumber, NFE> round = degree() == 2 ? new MathMap<>() {
					@Override
					public NFE evaluate(ComplexNumber t) {
						IntE alphaPart = r.divide(t.complexPart(), alpha.complexPart()).round();
						ComplexNumber residue = c.subtract(t, c.multiply(alphaPart, alpha));
						IntE realPart = residue.realPart().round();
						return add(getInteger(realPart), multiply(alphaPart, alpha()));
					}
				} : null;
				complexEmbeddings.add(new EmbeddedNumberField<>(this, c, alpha, round));
			}
		}
		return Collections.unmodifiableList(complexEmbeddings);
	}

	public List<EmbeddedNumberField<Ext, CompletedNumberField>> padicEmbeddings(PAdicField base) {
		IntE prime = base.getPrime();
		if (padicEmbeddings == null) {
			padicEmbeddings = new TreeMap<>();
		}
		if (!padicEmbeddings.containsKey(prime)) {
			padicEmbeddings.put(prime, new TreeMap<>());
		}
		int accuracy = base.getAccuracy();
		if (!padicEmbeddings.get(prime).containsKey(accuracy)) {
			List<EmbeddedNumberField<Ext, CompletedNumberField>> result = new ArrayList<>();
			for (NumberFieldIdeal primeIdeal : maximalOrder().idealsOver(prime)) {
				LocalizedNumberField localized = maximalOrder().localizeAndQuotient(primeIdeal);
				OtherVersion<NFE, Ext, FFE, CompletedNumberField> complete = localized.complete(accuracy);
				result.add(new EmbeddedNumberField<>(this, complete.getField(),
						complete.getRetraction().evaluate(alpha()), complete.getEmbedding()));
			}
//			UnivariatePolynomialRing<PAdicNumber> padicPolynomials = base.getUnivariatePolynomialRing();
//			Integers z = Integers.z();
//			DiscreteValuationRing<Fraction, PFE> localized = z.localize(prime);
//			PrimeField fp = PrimeField.getPrimeField(prime);
//			Extension<PFE, PFE, FFE, FiniteField> trivialExtension = fp
//					.getExtension(fp.getUnivariatePolynomialRing().getVar());
//			TheMontesResult<Fraction, PFE, PFE, FFE, FiniteField> theMontes = localized
//					.theMontesAlgorithm(minimalPolynomial(), trivialExtension);
//			for (OkutsuType<Fraction, PFE, PFE, FFE, FiniteField> type : theMontes.getTypes()) {
//				type = localized.singleFactorLifting(type, accuracy);
//				CompleteDVRExtension<PAdicNumber, PFE, PFE, FFE, FiniteField> embeddingField = CompleteDVRExtension
//						.getCompleteDVRExtension(
//								padicPolynomials.getEmbedding(type.representative(), base.fromRationalMap()), base,
//								trivialExtension);
//				Ext<PAdicNumber> alpha = embeddingField.alpha();
//				result.add(new EmbeddedNumberField<>(this, embeddingField, alpha, new MathMap<>() {
//					@Override
//					public NFE evaluate(Ext<PAdicNumber> t) {
//						t = embeddingField.round(t, accuracy);
//						UnivariatePolynomial<Fraction> asPolynomial = q.getUnivariatePolynomialRing()
//								.getEmbedding(t.asPolynomial(), base.toRationalMap());
//						UnivariatePolynomial<NFE> asNFEPolynomial = getUnivariatePolynomialRing()
//								.getEmbedding(asPolynomial, getEmbeddingMap());
//						return getUnivariatePolynomialRing().evaluate(asNFEPolynomial, alpha());
//					}
//				}));
//			}
			padicEmbeddings.get(prime).put(accuracy, result);
		}
		return padicEmbeddings.get(prime).get(accuracy);
	}

	public Real minkowskiBound() {
		Real discriminantSqrt = r.positiveSqrt(r.abs(r.getInteger(discriminant())));
		int complexEmbeddings = complexEmbeddings().size();
		Real fourOverPiToPower = r.power(r.divide(r.getInteger(4), r.pi()), complexEmbeddings);
		int degree = degree();
		Real nFactorialOverNToN = r.one();
		for (int i = 1; i <= degree; i++) {
			nFactorialOverNToN = r.multiply(nFactorialOverNToN, r.divide(r.getInteger(i), r.getInteger(degree)));
		}
		return r.multiply(discriminantSqrt, fourOverPiToPower, nFactorialOverNToN);
	}

	public Vector<Real> minkowskiEmbedding(NFE t) {
		List<Real> embedding = new ArrayList<>();
		for (EmbeddedNumberField<Real, Reals> realEmbedding : realEmbeddings()) {
			embedding.add(realEmbedding.embedding(t));
		}
		for (EmbeddedNumberField<ComplexNumber, Complex> complexEmbedding : complexEmbeddings()) {
			ComplexNumber embedded = complexEmbedding.embedding(t);
			embedding.add(embedded.realPart());
			embedding.add(embedded.complexPart());
		}
		return new Vector<>(embedding);
	}

	public FiniteRealVectorSpace minkowskiEmbeddingSpace() {
		if (minkowskiEmbeddingSpace == null) {
			minkowskiEmbeddingSpace = new FiniteRealVectorSpace(r, degree());
		}
		return minkowskiEmbeddingSpace;
	}

	public FiniteRealVectorSpace logRepresentationSpace() {
		if (logRepresentationSpace == null) {
			logRepresentationSpace = new FiniteRealVectorSpace(r, realEmbeddings().size() + complexEmbeddings().size());
		}
		return logRepresentationSpace;
	}

	public Vector<Real> logRepresentation(NFE t) {
		List<Real> result = new ArrayList<>();
		for (EmbeddedNumberField<Real, Reals> realEmbedding : realEmbeddings()) {
			result.add(r.log(realEmbedding.value(t)));
		}
		for (EmbeddedNumberField<ComplexNumber, Complex> complexEmbedding : complexEmbeddings()) {
			result.add(r.multiply(2, r.log(complexEmbedding.value(t))));
		}
		return new Vector<>(result);
	}

	public ExtendedResultantResult<NFE> extendedResultant(UnivariatePolynomial<NFE> t1, UnivariatePolynomial<NFE> t2) {
		UnivariatePolynomialRing<NFE> polynomials = getUnivariatePolynomialRing();
		if (t2.degree() == 0) {
			if (t1.degree() == 0) {
				Pair<NFE, NFE> bezout = bezoutIdentity(t1.leadingCoefficient(), t2.leadingCoefficient());
				return new ExtendedResultantResult<>(one(), polynomials.getEmbedding(bezout.getFirst()),
						polynomials.getEmbedding(bezout.getSecond()), polynomials.one(),
						polynomials.getEmbedding(bezout.getFirst()), polynomials.getEmbedding(bezout.getSecond()));
			}
			NFE power = power(t2.leadingCoefficient(), t1.degree() - 1);
			NFE resultant = multiply(power, t2.leadingCoefficient());
			ExtendedResultantResult<NFE> result = new ExtendedResultantResult<>(resultant, polynomials.zero(),
					polynomials.getEmbedding(power), polynomials.getEmbedding(resultant), polynomials.zero(),
					polynomials.getEmbedding(power));
			return result;
		}
		NFE unit1 = t1.leadingCoefficient();
		t1 = polynomials.normalize(t1);
		NFE unit2 = t2.leadingCoefficient();
		t2 = polynomials.normalize(t2);
		Integers z = Integers.z();
		IntE max = z.one();
		for (int i = 0; i <= Math.max(t1.degree(), t2.degree()); i++) {
			Fraction norm = norm(t1.univariateCoefficient(i));
			IntE limit = z.archimedeanValue(z.multiply(norm.getNumerator(), norm.getDenominator()));
			if (limit.compareTo(max) > 0) {
				max = limit;
			}
			norm = norm(t2.univariateCoefficient(i));
			limit = z.archimedeanValue(z.multiply(norm.getNumerator(), norm.getDenominator()));
			if (limit.compareTo(max) > 0) {
				max = limit;
			}
		}
		max = z.multiply(2, max);
		NumberFieldIntegers order = maximalOrder();
		Map<NumberFieldIdeal, ExtendedResultantResult<FFE>> perPrime = new TreeMap<>();
		int minGcdDegree = Math.max(t1.degree(), t2.degree());
		IntE currentNorm = z.one();
		for (IntE prime : z.setOfPrimes()) {
			if (currentNorm.compareTo(max) > 0) {
				break;
			}
			for (NumberFieldIdeal ideal : order.idealsOver(prime)) {
				LocalizedNumberField localized = order.localizeAndQuotient(ideal);
				if (localized.ramificationIndex() != 1) {
					continue;
				}
				if (!localized.ringOfIntegers().valuationOfUnivariatePolynomial(t1).equals(Value.ZERO)
						|| !localized.ringOfIntegers().valuationOfUnivariatePolynomial(t1).equals(Value.ZERO)) {
					continue;
				}
				UnivariatePolynomial<FFE> reduced1 = localized.ringOfIntegers().reduceUnivariatePolynomial(t1);
				UnivariatePolynomial<FFE> reduced2 = localized.ringOfIntegers().reduceUnivariatePolynomial(t2);
				ExtendedResultantResult<FFE> resultant = localized.ringOfIntegers().reduction()
						.getUnivariatePolynomialRing().extendedResultant(reduced1, reduced2);
				if (resultant.getGcd().degree() < minGcdDegree) {
					minGcdDegree = resultant.getGcd().degree();
					currentNorm = z.one();
					perPrime.clear();
				}
				if (resultant.getGcd().degree() > minGcdDegree) {
					continue;
				}
				currentNorm = z.multiply(prime, currentNorm);
				perPrime.put(ideal, resultant);
			}
		}
		List<NumberFieldIdeal> ideals = new ArrayList<>();
		ideals.addAll(perPrime.keySet());
		ChineseRemainderPreparation<NFE> preparation = order.prepareChineseRemainderTheorem(ideals);
		List<NFE> resultantMod = new ArrayList<>();
		for (NumberFieldIdeal ideal : ideals) {
			LocalizedNumberField localized = order.localizeAndQuotient(ideal);
			resultantMod.add(localized.ringOfIntegers().lift(perPrime.get(ideal).getResultant()));
		}
		NFE resultant = order.chineseRemainderTheorem(resultantMod, preparation);
		UnivariatePolynomial<NFE> resultantCoeff1;
		UnivariatePolynomial<NFE> resultantCoeff2;
		List<NFE> gcdCoeffs = new ArrayList<>();
		for (int i = 0; i <= minGcdDegree; i++) {
			List<NFE> mods = new ArrayList<>();
			for (NumberFieldIdeal ideal : ideals) {
				LocalizedNumberField localized = order.localizeAndQuotient(ideal);
				mods.add(localized.ringOfIntegers().lift(perPrime.get(ideal).getGcd().univariateCoefficient(i)));
			}
		}
		UnivariatePolynomial<NFE> gcd;
		UnivariatePolynomial<NFE> coeff1;
		UnivariatePolynomial<NFE> coeff2;
		return new ExtendedResultantResult<>(resultant, null, null, null, null, null);
	}

	@Override
	public FactorizationResult<Polynomial<NFE>, NFE> factorization(UnivariatePolynomial<NFE> t) {
		if (t.degree() == 0) {
			return new FactorizationResult<>(t.univariateCoefficient(0), Collections.emptySortedMap());
		}
		UnivariatePolynomialRing<NFE> ring = getUnivariatePolynomialRing();
		// extendedResultant(t, ring.derivative(t));
		if (t.degree() == 1) {
			NFE unit = t.leadingCoefficient();
			t = ring.normalize(t);
			return new FactorizationResult<>(unit, SingletonSortedMap.map(t, 1));
		}
		SortedMap<Polynomial<NFE>, Integer> result = new TreeMap<>();
		if (!hasIntegerMinimalPolynomial) {
			UnivariatePolynomialRing<NFE> integerRing = withIntegerMinimalPolynomial().getUnivariatePolynomialRing();
			UnivariatePolynomial<NFE> integer = integerRing.getEmbedding(t, toIntegerMinimalPolynomialMap());
			FactorizationResult<Polynomial<NFE>, NFE> integerResult = withIntegerMinimalPolynomial()
					.factorization(integer);
			for (Polynomial<NFE> factor : integerResult.primeFactors()) {
				result.put(ring.getEmbedding(factor, fromIntegerMinimalPolynomialMap()),
						integerResult.multiplicity(factor));
			}
			return new FactorizationResult<>(fromIntegerMinimalPolynomial(integerResult.getUnit()), result);
		}
		NFE unit = t.leadingCoefficient();
		t = removeDenominators(t);
		Iterator<IntE> primes = Integers.z().primes();
		for (int i = 0; i < 10; i++) {
			primes.next();
		}
		int count = 0;
		UnivariatePolynomial<NFE> derivative = ring.derivative(t);
		boolean isSquareFree = false;
		while (count < 5 && !isSquareFree) {
			IntE prime = primes.next();
			for (NumberFieldIdeal ideal : maximalOrder().idealsOver(prime)) {
				LocalizedNumberField localized = maximalOrder().localizeAndQuotient(ideal);
				if (localized.ramificationIndex() != 1) {
					continue;
				}
				if (!localized.ringOfIntegers().valuationOfUnivariatePolynomial(t).equals(Value.ZERO)
						|| !localized.ringOfIntegers().valuationOfUnivariatePolynomial(derivative).equals(Value.ZERO)) {
					continue;
				}
				UnivariatePolynomial<FFE> reduced = localized.ringOfIntegers().reduceUnivariatePolynomial(t);
				UnivariatePolynomial<FFE> reducedDerivative = localized.ringOfIntegers()
						.reduceUnivariatePolynomial(derivative);
				ExtendedResultantResult<FFE> resultant = localized.ringOfIntegers().reduction()
						.getUnivariatePolynomialRing().extendedResultant(reduced, reducedDerivative);
				if (resultant.getGcd().degree() == 0) {
					isSquareFree = true;
					break;
				}
				count++;
			}
		}
		FactorizationResult<Polynomial<NFE>, NFE> squareFree;
		if (isSquareFree) {
			squareFree = new FactorizationResult<>(one(), SingletonSortedMap.map(t, 1));
		} else {
			squareFree = ring.squareFreeFactorization(t);
		}
		for (Polynomial<NFE> squareFreeFactor : squareFree.primeFactors()) {
			for (Polynomial<NFE> factor : factorizeSquareFree(removeDenominators(squareFreeFactor))) {
				result.put(ring.normalize(factor), squareFree.multiplicity(squareFreeFactor));
			}
		}
		return new FactorizationResult<>(unit, result);
	}

	private UnivariatePolynomial<NFE> removeDenominators(Polynomial<NFE> t) {
		UnivariatePolynomialRing<NFE> ring = getUnivariatePolynomialRing();
		UnivariatePolynomial<NFE> p = ring.normalize(t);
		Integers z = Integers.z();
		IntE lcm = z.one();
		for (int i = 0; i <= t.degree(); i++) {
			Vector<Fraction> c = matrixAlgebra().multiply(maximalOrder().toIntegralBasisBaseChange(),
					asVector(p.univariateCoefficient(i)));
			for (int j = 0; j < degree(); j++) {
				lcm = z.lcm(lcm, c.get(j + 1).getDenominator());
			}
		}
		return ring.multiply(getInteger(lcm), p);
	}

	private class CombinedFactors {
		private SortedSet<Integer> usedIndeces;
		private Polynomial<Ext> combined;

		public CombinedFactors(SortedSet<Integer> usedIndeces, Polynomial<Ext> combined) {
			this.usedIndeces = usedIndeces;
			this.combined = combined;
		}
	}

	private List<Polynomial<NFE>> factorizeSquareFree(UnivariatePolynomial<NFE> t) {
		if (t.degree() == 0) {
			return Collections.singletonList(t);
		}
		if (t.degree() < 0) {
			throw new ArithmeticException("Cannot factor 0");
		}
		if (polynomialFactorizationCache.containsKey(t)) {
			return polynomialFactorizationCache.get(t);
		}
		Integers z = Integers.z();
		UnivariatePolynomial<NFE> f = t;
		BigInteger limit = BigInteger.ZERO;
		for (Monomial m : t.monomials()) {
			NFE c = t.coefficient(m);
			BigInteger l = norm(c).asInteger().getValue().abs();
			if (l.compareTo(limit) > 0) {
				limit = l;
			}
		}
		limit = limit.add(BigInteger.ONE);
		limit = limit.shiftLeft(1);
		if (idealForFactorization == null) {
			IntE denominators = z.one();
			for (NFE basis : maximalOrder().getModuleGenerators()) {
				for (int i = 0; i < degree(); i++) {
					denominators = z.lcm(basis.asPolynomial().univariateCoefficient(i).getDenominator(), denominators);
				}
			}
			IntE discriminant = discriminant();
			int degreeRatio = degree();
			IntE max = z.getInteger(discriminant.getValue().abs());
			if (max.compareTo(z.getInteger(20)) > 0) {
				max = z.getInteger(20);
			}
			for (IntE prime : z.setOfPrimes()) {
				if (prime.compareTo(max) > 0 && idealForFactorization != null) {
					break;
				}
				if (z.isDivisible(denominators, prime)) {
					continue;
				}
				if (z.isDivisible(discriminant, prime)) {
					continue;
				}
				for (NumberFieldIdeal ideal : maximalOrder().idealsOver(prime)) {
					int degreeRatioThisIdeal = degree()
							/ (ideal.type().ramificationIndex() * ideal.type().residueDegree());
					if (degreeRatioThisIdeal <= degreeRatio) {
						idealForFactorization = ideal;
						degreeRatio = degreeRatioThisIdeal;
					}
					if (degreeRatio == 1) {
						break;
					}
				}
				if (degreeRatio == 1) {
					break;
				}
			}
		}
		IntE prime = idealForFactorization.intGenerator();
		int requiredAccuracy = (degree()
				/ (idealForFactorization.type().ramificationIndex() * idealForFactorization.type().residueDegree()))
				* (int) (Math.log(limit.doubleValue()) / Math.log(prime.getValue().doubleValue())) + 2;
		EmbeddedNumberField<Ext, CompletedNumberField> embeddedNumberFieldHigh = null;
		for (EmbeddedNumberField<Ext, CompletedNumberField> embedded : padicEmbeddings(
				new PAdicField(prime, requiredAccuracy * t.degree()))) {
			if (embedded.embeddingField().degree() == idealForFactorization.type().ramificationIndex()
					* idealForFactorization.type().residueDegree()) {
				embeddedNumberFieldHigh = embedded;
				break;
			}
		}
		DiscreteValuationRing<Ext, FFE> zp = embeddedNumberFieldHigh.embeddingField().ringOfIntegers();
		UnivariatePolynomialRing<Ext> zpr = zp.getUnivariatePolynomialRing();
		UnivariatePolynomial<Ext> fp = zpr.getEmbedding(t, embeddedNumberFieldHigh.embeddingMap());
		List<Polynomial<Ext>> liftedFactors = new ArrayList<>();
		FactorizationResult<Polynomial<Ext>, Ext> padicFactorizationOkutsu = zp.factorization(zpr.normalize(fp));
		EmbeddedNumberField<Ext, CompletedNumberField> embeddedNumberField = embeddedNumberFieldHigh;
		for (Polynomial<Ext> factor : padicFactorizationOkutsu.primeFactors()) {
			if (factor.degree() > 0) {
				liftedFactors.add(embeddedNumberField.embeddingField().getUnivariatePolynomialRing()
						.getEmbedding(factor, new MathMap<>() {
							@Override
							public Ext evaluate(Ext t) {
								return embeddedNumberField.embeddingField().round(t, requiredAccuracy);
							}
						}));
			}
		}

		Map<Integer, Polynomial<Ext>> padicFactors = new TreeMap<>();
		for (int i = 0; i < liftedFactors.size(); i++) {
			padicFactors.put(i, liftedFactors.get(i));
		}
		List<CombinedFactors> padicCombinedFactors = new ArrayList<>();
		padicCombinedFactors.add(new CombinedFactors(Collections.emptySortedSet(), zpr.one()));
		List<Polynomial<NFE>> intFactors = new ArrayList<>();
		while (padicFactors.size() > 0) {
			CheckCombinationsResult result = checkCombinations(t, padicFactors, padicCombinedFactors,
					embeddedNumberField, requiredAccuracy);
			t = result.cofactor;
			intFactors.addAll(result.factors);
			for (int k : result.usedFactors) {
				padicFactors.remove(k);
			}
			padicCombinedFactors = result.combined;
		}
		if (t.degree() != 0) {
			throw new ArithmeticException("Factorization recombination failed!");
		}
		polynomialFactorizationCache.put(f, intFactors);
		return intFactors;
	}

	private static class CheckCombinationsResult {
		private List<Polynomial<NFE>> factors = new ArrayList<>();
		private UnivariatePolynomial<NFE> cofactor;
		private Set<Integer> usedFactors = new TreeSet<>();
		private List<CombinedFactors> combined = new ArrayList<>();
	}

	private CheckCombinationsResult checkCombinations(UnivariatePolynomial<NFE> t,
			Map<Integer, Polynomial<Ext>> padicFactors, List<CombinedFactors> padicCombinedFactors,
			EmbeddedNumberField<Ext, CompletedNumberField> qp, int accuracy) {
		CheckCombinationsResult result = new CheckCombinationsResult();
		UnivariatePolynomialRing<Ext> ring = qp.embeddingField().getUnivariatePolynomialRing();
		result.cofactor = t;
		for (int i : padicFactors.keySet()) {
			Polynomial<Ext> padicFactor = padicFactors.get(i);
			for (CombinedFactors padicCombinedFactor : padicCombinedFactors) {
				if (padicCombinedFactor.usedIndeces.size() != 0 && padicCombinedFactor.usedIndeces.first() <= i) {
					continue;
				}
				SortedSet<Integer> indeces = new TreeSet<>();
				indeces.addAll(padicCombinedFactor.usedIndeces);
				indeces.add(i);
				Polynomial<Ext> newCombined = ring.multiply(padicFactor, padicCombinedFactor.combined);
				CheckFactorResult cfr = checkFactor(result.cofactor, newCombined, qp, accuracy);
				if (cfr.foundFactor) {
					result.factors.add(cfr.factor);
					result.cofactor = cfr.cofactor;
					result.usedFactors.addAll(indeces);
					break;
				} else {
					result.combined.add(new CombinedFactors(indeces, newCombined));
				}
			}
		}
		return result;
	}

	private static class CheckFactorResult {
		private boolean foundFactor = false;
		private Polynomial<NFE> factor = null;
		private UnivariatePolynomial<NFE> cofactor = null;
	}

	private CheckFactorResult checkFactor(UnivariatePolynomial<NFE> t, Polynomial<Ext> potentialFactor,
			EmbeddedNumberField<Ext, CompletedNumberField> qp, int accuracy) {
		CheckFactorResult result = new CheckFactorResult();
		UnivariatePolynomialRing<NFE> ring = getUnivariatePolynomialRing();
		UnivariatePolynomialRing<Ext> qpRing = qp.embeddingField().getUnivariatePolynomialRing();
		potentialFactor = qpRing.multiply(
				qp.embeddingField().divide(qp.embedding(t.leadingCoefficient()), potentialFactor.leadingCoefficient()),
				potentialFactor);
		Polynomial<NFE> factor = ring.getEmbedding(potentialFactor, new MathMap<>() {
			@Override
			public NFE evaluate(Ext number) {
				return qp.embeddingField().roundToInteger(number, accuracy);
			}
		});
		factor = ring.contentFree(factor);
		QuotientAndRemainderResult<Polynomial<NFE>> qr = ring.quotientAndRemainder(t, factor);
		if (qr.getRemainder().equals(ring.zero())) {
			result.foundFactor = true;
			result.factor = factor;
			result.cofactor = ring.toUnivariate(qr.getQuotient());
		}
		return result;
	}

//	private List<Polynomial<NFE>> factorizeSquareFree(UnivariatePolynomial<NFE> t) {
//		UnivariatePolynomialRing<NFE> ring = getUnivariatePolynomialRing();
//		if (t.degree() == 1) {
//			return Collections.singletonList(ring.normalize(t));
//		}
//		Rationals q = Rationals.q();
//		Integers z = Integers.z();
//		UnivariatePolynomialRing<Fraction> rationalRing = q.getUnivariatePolynomialRing();
//		GenericAlgebraicRingExtension<NFE> cr = new GenericAlgebraicRingExtension<>(t, this);
//		Iterator<UnivariatePolynomial<IntE>> it = z.getUnivariatePolynomialRing().polynomials(degree() - 1);
//		while (true) {
//			NFE probeTerm = fromSmallDegreePolynomial(rationalRing.getEmbedding(it.next(), new MathMap<>() {
//				@Override
//				public Fraction evaluate(IntE t) {
//					return q.getEmbedding(t);
//				}
//			}));
//			GenericAlgebraicExtensionElement<NFE> probe = cr.add(cr.alpha(),
//					cr.getEmbedding(multiply(probeTerm, alpha())));
//			List<Vector<Fraction>> generators = new ArrayList<>();
//			GenericAlgebraicExtensionElement<NFE> power = cr.one();
//			for (int i = 0; i < t.degree() * degree(); i++) {
//				generators.add(asRationalVector(power, t.degree()));
//				power = cr.multiply(power, probe);
//			}
//			FiniteVectorSpace<Fraction> asVectorSpace = new FiniteVectorSpace<>(q, t.degree() * degree());
//			Matrix<Fraction> m = Matrix.fromColumns(generators);
//			MatrixAlgebra<Fraction> algebra = asVectorSpace.matrixAlgebra();
//			if (algebra.rank(m) != t.degree() * degree()) {
//				continue;
//			}
//			UnivariatePolynomial<Fraction> x = rationalRing.getPolynomial(
//					algebra.solve(m, asRationalVector(cr.fromPolynomial(ring.getVar()), t.degree())).asList());
//			UnivariatePolynomial<Fraction> gamma = rationalRing
//					.getPolynomial(algebra.solve(m, asRationalVector(cr.getEmbedding(alpha()), t.degree())).asList());
//			Vector<Fraction> lastPower = asRationalVector(power, t.degree());
//			UnivariatePolynomial<Fraction> minimalPolynomial = rationalRing
//					.getPolynomial(algebra.solve(m, lastPower).asList());
//			minimalPolynomial = rationalRing.toUnivariate(
//					rationalRing.subtract(rationalRing.getVarPower(t.degree() * degree()), minimalPolynomial));
//			FactorizationResult<Polynomial<Fraction>, Fraction> rationalFactors = q.factorization(minimalPolynomial);
//			if (rationalFactors.primeFactors().size() == 1) {
//				return Collections.singletonList(ring.normalize(t));
//			}
//			List<Polynomial<NFE>> factors = new ArrayList<>();
//			for (Polynomial<Fraction> factor : rationalFactors.primeFactors()) {
//				NumberField extension = new NumberField(rationalRing.toUnivariate(factor));
//				if (!extension.minimalPolynomial().equals(factor)) {
//					throw new ArithmeticException("Something messed with the minimalpolynomial!");
//				}
//				UnivariatePolynomialRing<NFE> extensionRing = extension.getUnivariatePolynomialRing();
//				NFE embeddedX = extensionRing.evaluate(extensionRing.getEmbedding(x, extension.getEmbeddingMap()),
//						extension.alpha());
//				NFE embeddedGamma = extensionRing
//						.evaluate(extensionRing.getEmbedding(gamma, extension.getEmbeddingMap()), extension.alpha());
//				List<Vector<Fraction>> gammaXBase = new ArrayList<>();
//				NFE xPower = extension.one();
//				for (int i = 0; i < extension.degree() / degree(); i++) {
//					NFE gammaPower = extension.one();
//					for (int j = 0; j < degree(); j++) {
//						gammaXBase.add(extension.asVector(extension.multiply(xPower, gammaPower)));
//						gammaPower = extension.multiply(gammaPower, embeddedGamma);
//					}
//					xPower = extension.multiply(xPower, embeddedX);
//				}
//				Matrix<Fraction> gammaXMatrix = Matrix.fromColumns(gammaXBase);
//				UnivariatePolynomial<NFE> nfeFactor = minimalPolynomial(embeddedX, extension.degree() / degree(),
//						extension, this, new MathMap<>() {
//
//							@Override
//							public Vector<NFE> evaluate(NFE t) {
//								Vector<Fraction> overGammaX = extension.matrixAlgebra().solve(gammaXMatrix,
//										extension.asVector(t));
//								List<NFE> asList = new ArrayList<>();
//								for (int i = 0; i < extension.degree() / degree(); i++) {
//									asList.add(fromPolynomial(rationalRing.getPolynomial(
//											overGammaX.asList().subList(i * degree(), (i + 1) * degree()))));
//								}
//								return new Vector<>(asList);
//							}
//						});
//				factors.add(ring.normalize(nfeFactor));
//			}
//			return factors;
//		}
//	}

	private Vector<Fraction> asRationalVector(GenericAlgebraicExtensionElement<NFE> e, int dimension) {
		UnivariatePolynomialRing<NFE> ring = getUnivariatePolynomialRing();
		List<NFE> asNFEList = ring.asVector(e.asPolynomial(), dimension - 1).asList();
		List<Fraction> asList = new ArrayList<>();
		for (NFE nfe : asNFEList) {
			asList.addAll(asVector(nfe).asList());
		}
		return new Vector<>(asList);
	}

	public UnivariatePolynomial<NFE> integerMinimalPolynomial(Complex c, List<ComplexNumber> conjugates) {
		Integers z = Integers.z();
		if (this.degree() != 2 || this.discriminant().compareTo(z.zero()) > 0) {
			throw new ArithmeticException(
					"Currently only implemented for imaginary quadratic number fields because it's fiddly otherwise!");
		}
		Reals r = c.getReals();
		int degree = conjugates.size();
		ComplexNumber complexGenerator = complexEmbeddings().get(0)
				.embedding(maximalOrder().getModuleGenerators().get(1));
		List<ComplexNumber> powerValues = new ArrayList<>();
		List<ComplexNumber> multipliedPowerValues = new ArrayList<>();
		for (int i = 0; i < degree; i++) {
			powerValues.add(c.one());
			multipliedPowerValues.add(complexGenerator);
		}
		List<Vector<Real>> powers = new ArrayList<>();
		for (int i = 0; i < degree; i++) {
			powers.add(asEmbeddingVector(powerValues));
			powers.add(asEmbeddingVector(multipliedPowerValues));
			for (int j = 0; j < degree; j++) {
				powerValues.set(j, c.multiply(conjugates.get(j), powerValues.get(j)));
				multipliedPowerValues.set(j, c.multiply(complexGenerator, powerValues.get(j)));
			}
		}
		Vector<Real> target = asEmbeddingVector(powerValues);
		FiniteRealVectorSpace space = new FiniteRealVectorSpace(r, 2 * degree);
		RealLattice lattice = new RealLattice(space, powers);
		Vector<IntE> integerTarget = lattice.asVector(target);
		List<Vector<IntE>> integerPowers = new ArrayList<>();
		for (int i = 0; i < 2 * degree; i++) {
			integerPowers.add(lattice.asVector(powers.get(i)));
		}
		Matrix<IntE> integerBasis = Matrix.fromColumns(integerPowers);
		Vector<IntE> solvedIntegerVector = integerBasis.getModule(z).solve(integerBasis, integerTarget);
		List<NFE> coefficients = new ArrayList<>();
		for (int i = 0; i < degree; i++) {
			List<IntE> coeff = new ArrayList<>();
			coeff.add(solvedIntegerVector.get(2 * i + 1));
			coeff.add(solvedIntegerVector.get(2 * i + 2));
			coefficients.add(negative(maximalOrder().fromVector(new Vector<>(coeff))));
		}
		coefficients.add(one());
		return getUnivariatePolynomialRing().getPolynomial(coefficients);
	}

	private Vector<Real> asEmbeddingVector(List<ComplexNumber> conjugates) {
		List<Real> values = new ArrayList<>();
		for (ComplexNumber conjugate : conjugates) {
			values.add(conjugate.realPart());
			values.add(conjugate.complexPart());
		}
		return new Vector<>(values);
	}

	public Pair<UnivariatePolynomial<NFE>, UnivariatePolynomial<NFE>> simplifyMinimalPolynomial(
			UnivariatePolynomial<NFE> minimalPolynomial) {
		UnivariatePolynomialRing<NFE> polynomialRing = getUnivariatePolynomialRing();
		minimalPolynomial = polynomialRing.normalize(minimalPolynomial);
		boolean rational = true;
		for (int i = 0; i <= minimalPolynomial.degree(); i++) {
			rational = rational && minimalPolynomial.univariateCoefficient(i).asPolynomial().degree() <= 0;
		}
		if (minimalPolynomial.degree() == 1) {
			return new Pair<>(getUnivariatePolynomialRing().getVar(), getUnivariatePolynomialRing().getEmbedding(divide(
					negative(minimalPolynomial.univariateCoefficient(0)), minimalPolynomial.univariateCoefficient(1))));
		}
		if (rational) {
			UnivariatePolynomialRing<Fraction> rationalPolynomialRing = q.getUnivariatePolynomialRing();
			UnivariatePolynomial<Fraction> rationalMinimalPolynomial = rationalPolynomialRing
					.getEmbedding(minimalPolynomial, new MathMap<>() {
						@Override
						public Fraction evaluate(NFE t) {
							return asVector(t).get(1);
						}
					});
			NumberField nf = NumberField.getNumberField(rationalMinimalPolynomial);
			for (NFE t : nf.maximalOrder().getModuleGenerators()) {
				UnivariatePolynomial<Fraction> result = nf.minimalPolynomial(t);
				if (result.degree() == minimalPolynomial.degree()) {
					List<Vector<Fraction>> powers = new ArrayList<>();
					for (int i = 0; i < result.degree(); i++) {
						powers.add(nf.asVector(nf.power(t, i)));
					}
					Matrix<Fraction> baseChange = Matrix.fromColumns(powers);
					Matrix<Fraction> inverseBaseChange = nf.matrixAlgebra().inverse(baseChange);
					UnivariatePolynomial<Fraction> inSimplified = q.getUnivariatePolynomialRing()
							.getPolynomial(inverseBaseChange.column(2).asList());
					return new Pair<>(polynomialRing.getEmbedding(result, getEmbeddingMap()),
							polynomialRing.getEmbedding(inSimplified, getEmbeddingMap()));
				}
			}
		}
		Pair<UnivariatePolynomial<NFE>, UnivariatePolynomial<NFE>> denominatorSubstition = substituteDenominator(
				minimalPolynomial);
		minimalPolynomial = denominatorSubstition.getFirst();
		UnivariatePolynomial<NFE> inSimplified = denominatorSubstition.getSecond();
		NumberFieldIntegers order = maximalOrder();
		IdealGroup idealGroup = idealGroup();
		Iterator<Ideal<NFE>> it = null;
		NFE discriminant = null;
		ValueFractions q = Rationals.q().withInfValue();
		while (true) {
			if (it == null) {
				discriminant = polynomialRing.discriminant(minimalPolynomial);
				FractionalIdeal discrimantIdeal = idealGroup.getPrincipalIdeal(discriminant);
				FactorizationResult<Ideal<NFE>, Ideal<NFE>> numeratorFactorization = order
						.idealFactorization(discrimantIdeal.getNumerator());
				FactorizationResult<Ideal<NFE>, Ideal<NFE>> denominatorFactorization = order
						.idealFactorization(discrimantIdeal.getDenominator());
				SortedSet<Ideal<NFE>> primes = new TreeSet<>();
				primes.addAll(numeratorFactorization.primeFactors());
				primes.addAll(denominatorFactorization.primeFactors());
				it = primes.iterator();
			}
			if (!it.hasNext()) {
				break;
			}
			NumberFieldIdeal ideal = (NumberFieldIdeal) it.next();
			LocalizedNumberField localized = order.localizeAndQuotient(ideal);
			UnivariatePolynomialRing<FFE> reductionPolynomialRing = localized.residueField()
					.getUnivariatePolynomialRing();
			UnivariatePolynomial<FFE> reduced = localized.ringOfIntegers()
					.reduceUnivariatePolynomial(minimalPolynomial);
			FactorizationResult<Polynomial<FFE>, FFE> factorization = localized.residueField().factorization(reduced);
			UnivariatePolynomial<FFE> factor = reductionPolynomialRing.normalize(factorization.firstPrimeFactor());
			if (!factorization.isPrimePower() || factorization.firstPrimeFactor().degree() != 1) {
				continue;
			}
			NFE lift = localized.liftToInteger(factor.univariateCoefficient(0));
			UnivariatePolynomial<NFE> inverted = polynomialRing.getPolynomial(negative(lift), one());
			UnivariatePolynomial<NFE> newMinimalPolynomial = polynomialRing.substitute(minimalPolynomial,
					Collections.singletonList(inverted));
			UnivariatePolynomial<NFE> newInSimplified = polynomialRing.substitute(inSimplified,
					Collections.singletonList(polynomialRing.getPolynomial(lift, one())));
			Pair<UnivariatePolynomial<NFE>, UnivariatePolynomial<NFE>> integral = localized.ringOfIntegers()
					.integralPolynomial(newMinimalPolynomial, getInteger(ideal.intGenerator()));
			newMinimalPolynomial = integral.getFirst();
			newInSimplified = polynomialRing.substitute(newInSimplified,
					Collections.singletonList(integral.getSecond()));
			denominatorSubstition = substituteDenominator(newMinimalPolynomial);
			newMinimalPolynomial = denominatorSubstition.getFirst();
			newInSimplified = polynomialRing.substitute(newInSimplified,
					Collections.singletonList(denominatorSubstition.getSecond()));
			NFE newDiscriminant = polynomialRing.discriminant(newMinimalPolynomial);
			if (q.value(norm(newDiscriminant)).compareTo(q.value(norm(discriminant))) < 0) {
				minimalPolynomial = newMinimalPolynomial;
				inSimplified = newInSimplified;
				it = null;
			}
		}
		return new Pair<>(polynomialRing.depress(polynomialRing.normalize(minimalPolynomial)), null);
	}

	private Pair<UnivariatePolynomial<NFE>, UnivariatePolynomial<NFE>> substituteDenominator(
			UnivariatePolynomial<NFE> t) {
		UnivariatePolynomialRing<NFE> polynomialRing = getUnivariatePolynomialRing();
		t = polynomialRing.normalize(t);
		Integers z = Integers.z();
		IntE denominator = z.one();
		for (int i = 0; i < t.degree(); i++) {
			IntE denom = z.one();
			for (int j = 0; j < degree(); j++) {
				denom = z.lcm(t.univariateCoefficient(i).asPolynomial().univariateCoefficient(j).getDenominator(),
						denom);
			}
			FactorizationResult<IntE, IntE> factors = z.uniqueFactorization(denom);
			for (IntE prime : factors.primeFactors()) {
				denominator = z.lcm(
						z.power(prime, MiscAlgorithms.DivRoundUp(factors.multiplicity(prime), t.degree() - i)),
						denominator);
			}
		}
		UnivariatePolynomial<NFE> undo = polynomialRing.getEmbedding(getInteger(denominator), 1);
		t = polynomialRing.substitute(t,
				Collections.singletonList(polynomialRing.getEmbedding(inverse(getInteger(denominator)), 1)));
		t = polynomialRing.normalize(t);
		return new Pair<>(t, undo);
	}

	public IntE discriminant() {
		return maximalOrder().discriminant(maximalOrder().getModuleGenerators());
	}

	public BigInteger classNumber() {
		return idealClassGroup().getNumberOfElements();
	}

	public IdealGroup idealGroup() {
		if (idealGroup == null) {
			idealGroup = new IdealGroup(maximalOrder());
		}
		return idealGroup;
	}

	public IdealClassGroup idealClassGroup() {
		if (idealClassGroup == null) {
			idealClassGroup = new IdealClassGroup(maximalOrder());
		}
		return idealClassGroup;
	}

	public FieldEmbedding<Fraction, NFE, NumberField> hilbertClassField() {
		Integers z = Integers.z();
		if (degree() == 1) {
			return getEmbeddedExtension(getUnivariatePolynomialRing().getVar());
		}
		if (degree() == 2 && discriminant().compareTo(z.zero()) < 0) {
			return EllipticCurve.computeHilbertClassField(this).getClassField();
		}
		IntE classNumber = z.getInteger(classNumber());
		System.out.println("Class number: " + classNumber);
		FactorizationResult<IntE, IntE> factorization = z.uniqueFactorization(classNumber);
		List<UnivariatePolynomial<NFE>> resultMinimalPolynomials = new ArrayList<>();
		for (IntE prime : factorization.primeFactors()) {
			System.out.println("Finding sub field of order " + prime);
			// TODO
			return hilbertClassField(prime);
		}
		throw new ArithmeticException("Not implemented");
	}

	public FieldEmbedding<Fraction, NFE, NumberField> hilbertClassField(IntE prime) {
		System.out.println("Computing Hilbert Class Field of " + this + " (" + prime + ")");
		Integers z = Integers.z();
		UnivariatePolynomialRing<NFE> polynomialRing = getUnivariatePolynomialRing();
		IdealClassGroup icg = idealClassGroup();
		UnivariatePolynomial<NFE> unity = polynomialRing.toUnivariate(polynomialRing.divideChecked(
				polynomialRing.subtract(polynomialRing.getVarPower(prime.intValueExact()), polynomialRing.one()),
				polynomialRing.subtract(polynomialRing.getVar(), polynomialRing.one())));
		FactorizationResult<Polynomial<NFE>, NFE> rootsOfUnity = factorization(unity);
		System.out.println("Roots of unity " + (unity.degree() > 1 && rootsOfUnity.isIrreducible() ? "not " : "")
				+ "contained in field");
		if (unity.degree() > 1 && rootsOfUnity.isIrreducible()) {
			FieldEmbedding<Fraction, NFE, NumberField> embedding = getEmbeddedExtension(unity);
			return embedding.getField().hilbertClassField(prime);
		}
		NumberFieldIntegers order = maximalOrder();
		NumberFieldIdeal ideal = order.getIdeal(Collections.singletonList(getInteger(prime)));
		List<NFE> unitGroupGenerators = order.unitGroupGenerators();
		System.out.println("Unit Group Generators: " + unitGroupGenerators);
		List<IdealClass> idealClassGroupGenerators = icg.getModuleGenerators();
		System.out.println("Ideal Class Group Generators: " + idealClassGroupGenerators);
		List<Integer> ranks = idealClassGroup.getOrders();
		List<NFE> alpha = new ArrayList<>();
		for (int i = 0; i < idealClassGroupGenerators.size(); i++) {
			IdealClass ic = idealClassGroupGenerators.get(i);
			IntE rank = z.getInteger(ranks.get(i));
			if (!z.isDivisible(rank, prime)) {
				continue;
			}
			NFE multiplier = one();
			for (Ideal<NFE> primeIdeal : order.idealFactorization(ideal).primeFactors()) {
				int multiplicity = order.idealFactorization(ic.representative()).multiplicity(primeIdeal);
				if (multiplicity != 0) {
					multiplier = multiply(power(((NumberFieldIdeal) primeIdeal).uniformizer(), multiplicity),
							multiplier);
				}
			}
			NumberFieldIdeal power = order.power(ic.representative(), ranks.get(i));
			NFE generator = power.principalGenerator();
			alpha.add(divide(generator, power(multiplier, rank)));
		}
		FactorizationResult<Ideal<NFE>, Ideal<NFE>> idealFactorization = order.idealFactorization(ideal);
		System.out.println(ideal + " = " + idealFactorization);
		System.out.println("Resolved Ideal Classes: " + alpha);
		List<Vector<IntE>> columns = new ArrayList<>();
		for (Ideal<NFE> idealFactor : idealFactorization.primeFactors()) {
			NumberFieldIdeal primeIdealFactor = (NumberFieldIdeal) idealFactor;
			System.out.println("Processing: " + idealFactor);
			LocalizedNumberField localized = order.localizeAndQuotient(primeIdealFactor);
			ModuloNumberFieldIdeal mod = order
					.power(primeIdealFactor, prime.intValueExact() * idealFactorization.multiplicity(idealFactor))
					.modOut();
			System.out.println(mod + "* = <" + mod.getUnitGenerators() + ">");
			for (int i = 0; i < unitGroupGenerators.size(); i++) {
				NFE unitGenerator = unitGroupGenerators.get(i);
				Vector<IntE> column = mod.asUnitGeneratorVector(mod.reduce(unitGenerator));
				System.out.println("Unit: " + unitGenerator + " = " + column);
				if (columns.size() == i) {
					columns.add(column);
				} else {
					List<IntE> extendedColumn = new ArrayList<>();
					extendedColumn.addAll(columns.get(i).asList());
					extendedColumn.addAll(column.asList());
					columns.set(i, new Vector<>(extendedColumn));
				}
			}
			for (int i = 0; i < alpha.size(); i++) {
				NFE alphaGenerator = localized.divide(alpha.get(i),
						localized.power(localized.uniformizer(), localized.valuation(alpha.get(i)).value()));
				alphaGenerator = localized.round(alphaGenerator,
						prime.intValueExact() * idealFactorization.multiplicity(idealFactor));
				Vector<IntE> column = mod.asUnitGeneratorVector(mod.reduce(alphaGenerator));
				System.out.println("Alpha: " + alphaGenerator + " = " + column);
				int index = unitGroupGenerators.size() + i;
				if (columns.size() == index) {
					columns.add(column);
				} else {
					List<IntE> extendedColumn = new ArrayList<>();
					extendedColumn.addAll(columns.get(index).asList());
					extendedColumn.addAll(column.asList());
					columns.set(index, new Vector<>(extendedColumn));
				}
			}
		}
		PrimeField reduction = PrimeField.getPrimeField(prime);
		List<Vector<PFE>> reduced = new ArrayList<>();
		for (Vector<IntE> column : columns) {
			List<PFE> reducedColumn = new ArrayList<>();
			for (IntE value : column.asList()) {
				reducedColumn.add(reduction.reduce(value));
			}
			reduced.add(new Vector<>(reducedColumn));
		}
		List<NFE> generators = new ArrayList<>();
		generators.addAll(unitGroupGenerators);
		generators.addAll(alpha);
		Matrix<PFE> reducedMatrix = Matrix.fromColumns(reduced);
		System.out.println("Matrix:");
		System.out.println(reducedMatrix);
		List<Vector<PFE>> kernelBasis = reducedMatrix.getModule(reduction).kernelBasis(reducedMatrix);
		System.out.println("Kernel of matrix: " + kernelBasis);
		FieldEmbedding<Fraction, NFE, NumberField> result = getEmbeddedExtension(
				getUnivariatePolynomialRing().getVar());
		for (Vector<PFE> basisVector : kernelBasis) {
			NFE t = order.one();
			for (int i = 0; i < generators.size(); i++) {
				t = order.multiply(order.power(generators.get(i), reduction.liftToInteger(basisVector.get(i + 1))), t);
			}
			UnivariatePolynomialRing<NFE> resultPolynomialRing = result.getField().getUnivariatePolynomialRing();
			UnivariatePolynomial<NFE> kummer = resultPolynomialRing.toUnivariate(resultPolynomialRing.subtract(
					resultPolynomialRing.getVarPower(prime.intValueExact()), resultPolynomialRing.getEmbedding(t)));
			if (result.getField().isIrreducible(kummer)) {
				result = new FieldEmbedding<>(result, result.getField().getEmbeddedExtension(kummer));
			}
		}
		return result;
	}

	public MathMap<ComplexNumber, ComplexNumber> dirichletZetaFunction(int precision) {
		Real epsilon = r.power(r.getInteger(2), -precision);
		return new MathMap<>() {
			@Override
			public ComplexNumber evaluate(ComplexNumber t) {
				if (t.realPart().compareTo(r.one()) <= 0) {
					return null;
				}
				Integers z = Integers.z();
				Fraction realPartFraction = null;
				boolean useFraction = t.complexPart().equals(r.zero());
				if (useFraction) {
					realPartFraction = r.roundToFraction(t.realPart(), PRECISION_FOR_EMBEDDINGS);
					if (realPartFraction.getDenominator().compareTo(z.getInteger(10000)) > 0) {
						useFraction = false;
					}
				}
				ComplexNumber result = c.one();// .subtract(t, c.one());
				ComplexNumber prevResult;
				Iterator<IntE> primeIt = z.primes();
				do {
					prevResult = result;
					IntE prime = primeIt.next();
					System.out.println(result);
					System.out.println(prime);
					List<NumberFieldIdeal> ideals = maximalOrder().idealsOver(prime);
					for (NumberFieldIdeal ideal : ideals) {
						ComplexNumber multiplier;
						if (useFraction) {
							ComplexNumber coefficient = c.getEmbedding(
									q.power(q.getEmbedding(ideal.norm()), z.negative(realPartFraction.getNumerator())));
							int power = realPartFraction.denominatorIntValueExact();
							ComplexNumber secondCoefficient = c.multiply(power, coefficient);
							multiplier = c.getInteger(2);
							ComplexNumber prevMultiplier;
							do {
								prevMultiplier = multiplier;
								ComplexNumber eval = c.subtract(c.power(c.subtract(multiplier, c.one()), power),
										c.multiply(coefficient, c.power(multiplier, power)));
								ComplexNumber derivativeEval = c.subtract(
										c.multiply(power, c.power(c.subtract(multiplier, c.one()), power - 1)),
										c.multiply(secondCoefficient, c.power(multiplier, power - 1)));
								multiplier = c.subtract(multiplier, c.divide(eval, derivativeEval));
							} while (!c.close(multiplier, prevMultiplier));
						} else {
							ComplexNumber normToMinusS = c
									.exp(c.multiply(-1, c.getEmbedding(r.log(r.getInteger(ideal.norm()))), t));
							ComplexNumber denominator = c.subtract(c.one(), normToMinusS);
							multiplier = c.inverse(denominator);
						}
						result = c.multiply(result, multiplier);
					}
				} while (c.value(c.subtract(prevResult, result)).compareTo(epsilon) > 0);
				return result;
			}
		};
	}

	@Override
	public GaloisGroup<Fraction, NFE, NumberField> galoisGroup() {
		if (!hasIntegerMinimalPolynomial) {
			return withIntegerMinimalPolynomial().galoisGroup();
		}
		if (galoisGroup == null) {
			if (degree() > 7) {
				throw new ArithmeticException("No. Just no.");
			}
			Integers z = Integers.z();
			UnivariatePolynomialRing<IntE> univariate = z.getUnivariatePolynomialRing();
			PolynomialRing<IntE> multivariate = AbstractPolynomialRing.getPolynomialRing(z, degree(), Monomial.GREVLEX);
			List<Polynomial<IntE>> rootPolynomials = new ArrayList<>();
			if (degree() == 1) {
				System.out.println("C1");
				return null;
			}
			if (degree() == 2) {
				System.out.println("C2");
				return null;
			}
			Polynomial<IntE> discriminant = multivariate.one();
			for (int i = 1; i < degree(); i++) {
				for (int j = 0; j < i; j++) {
					discriminant = multivariate.multiply(
							multivariate.subtract(multivariate.getVar(i + 1), multivariate.getVar(j + 1)),
							discriminant);
				}
			}
			Set<Polynomial<IntE>> discriminantOrbit = new TreeSet<>();
			discriminantOrbit.add(discriminant);
			discriminantOrbit.add(multivariate.negative(discriminant));
			if (degree() == 4) {
				rootPolynomials.add(multivariate.add(multivariate.getVar(1), multivariate.getVar(2)));
				rootPolynomials.add(multivariate.subtract(multivariate.getVar(1), multivariate.getVar(2)));
			} else if (degree() == 5) {
				rootPolynomials.add(multivariate.subtract(multivariate.getVar(1), multivariate.getVar(2)));
				rootPolynomials.add(multivariate
						.power(multivariate.subtract(multivariate.add(multivariate.getVar(1), multivariate.getVar(2)),
								multivariate.add(multivariate.getVar(3), multivariate.getVar(4))), 2));

			} else if (degree() == 6) {
				rootPolynomials.add(multivariate.add(multivariate.getVar(1), multivariate.getVar(2)));
				rootPolynomials
						.add(multivariate.add(multivariate.getVar(1), multivariate.getVar(2), multivariate.getVar(3)));
				rootPolynomials.add(multivariate.subtract(multivariate.getVar(1), multivariate.getVar(2)));
				rootPolynomials.add(multivariate.add(multivariate.getVar(1), multivariate.getVar(2),
						multivariate.getVar(3), discriminant));
			} else if (degree() == 7) {
				rootPolynomials
						.add(multivariate.add(multivariate.getVar(1), multivariate.getVar(2), multivariate.getVar(3)));
			}
			UnivariatePolynomial<IntE> minimalPolynomial = univariate.getEmbedding(minimalPolynomial(),
					Rationals.q().getAsIntegerMap());
			UnivariatePolynomial<IntE> discriminantResolvant = z.tschirnhausenResolvant(minimalPolynomial,
					discriminantOrbit);
			System.out.println("Delta = " + z.factorization(discriminantResolvant));
			for (Polynomial<IntE> roots : rootPolynomials) {
				UnivariatePolynomial<IntE> resolvant = z.tschirnhausenResolvant(minimalPolynomial, roots);
				System.out.println(roots + " " + z.factorization(resolvant));
			}
//			if (isNormal()) {
//				List<FieldAutomorphism<Fraction, NFE, NumberField>> result = new ArrayList<>();
//				UnivariatePolynomialRing<Fraction> basePolynomials = q.getUnivariatePolynomialRing();
//				List<NFE> conjugates = conjugates(alpha());
//				for (NFE conjugate : conjugates) {
//					int[] map = new int[degree()];
//					for (int i = 0; i < conjugates.size(); i++) {
//						NFE image = fromPolynomial(basePolynomials.substitute(conjugates.get(i).asPolynomial(),
//								Collections.singletonList(conjugate.asPolynomial())));
//					map[i] =	conjugates.indexOf(image);
//					}
//					result.add(new FieldAutomorphism<>(this, map));
//				}
//				galoisGroup = new GaloisGroup<>(this, result);
//			}
		}
		return galoisGroup;
	}

	@Override
	protected NFE fromSmallDegreePolynomial(UnivariatePolynomial<Fraction> t) {
		return new NFE(this, t);
	}

	@Override
	protected NumberField asExtensionType() {
		return this;
	}

	public Fraction traceForm(NFE s1, NFE s2) {
		return trace(multiply(s1, s2));
	}

}
