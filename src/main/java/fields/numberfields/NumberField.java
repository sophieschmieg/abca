package fields.numberfields;

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.SortedMap;
import java.util.TreeMap;

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
import fields.helper.GenericAlgebraicRingExtension;
import fields.helper.GenericAlgebraicRingExtension.GenericAlgebraicExtensionElement;
import fields.integers.Integers;
import fields.integers.Integers.IntE;
import fields.integers.Rationals;
import fields.integers.Rationals.Fraction;
import fields.interfaces.AlgebraicExtensionElement;
import fields.interfaces.DiscreteValuationRing;
import fields.interfaces.DiscreteValuationRing.OkutsuType;
import fields.interfaces.DiscreteValuationRing.TheMontesResult;
import fields.interfaces.MathMap;
import fields.interfaces.Polynomial;
import fields.interfaces.UnivariatePolynomial;
import fields.interfaces.UnivariatePolynomialRing;
import fields.local.CompleteDVRExtension;
import fields.local.CompleteDVRExtension.Ext;
import fields.local.PAdicField;
import fields.local.PAdicField.PAdicNumber;
import fields.numberfields.IdealClassGroup.IdealClass;
import fields.numberfields.NumberField.NFE;
import fields.numberfields.NumberFieldIntegers.NumberFieldIdeal;
import fields.vectors.FiniteVectorSpace;
import fields.vectors.Matrix;
import fields.vectors.MatrixAlgebra;
import fields.vectors.Vector;
import util.MiscAlgorithms;

public class NumberField extends AbstractFieldExtension<Fraction, NFE, NumberField> {
	private Rationals q;
	private Reals r;
	private Complex c;
	private boolean embeddingsComputed;
	private List<EmbeddedNumberField<Real>> realEmbeddings;
	private List<EmbeddedNumberField<ComplexNumber>> complexEmbeddings;
	private FiniteRealVectorSpace logRepresentationSpace;
	private FiniteRealVectorSpace minkowskiEmbeddingSpace;
	private NumberFieldIntegers maximalOrder;
	private IdealGroup idealGroup;
	private IdealClassGroup idealClassGroup;
	private static final int PRECISION_FOR_EMBEDDINGS = 128;

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
			return this.e.toString("α", true);
		}

		@Override
		public UnivariatePolynomial<Fraction> asPolynomial() {
			return e;
		}

	}

	public NumberField() {
		super(Rationals.q());
		init();
	}

	public NumberField(UnivariatePolynomial<Fraction> minimalPolynomial) {
		super(minimalPolynomial, Rationals.q());
		init();
	}

	private void init() {
		this.q = Rationals.q();
		this.realEmbeddings = new ArrayList<>();
		this.complexEmbeddings = new ArrayList<>();
		this.embeddingsComputed = false;
		this.r = Reals.r(PRECISION_FOR_EMBEDDINGS);
		this.c = Complex.c(PRECISION_FOR_EMBEDDINGS);
		if (!q.isIrreducible(minimalPolynomial())) {
			throw new ArithmeticException("Not irreducible");
		}
	}

	public NumberFieldIntegers maximalOrder() {
		if (maximalOrder == null) {
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
		return new NumberField(minimalPolynomial);
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
			if (alpha.complexPart().equals(r.zero()) || alpha.complexPart().exponent() < -r.precision() + 5) {
				MathMap<Real, NFE> round = degree() == 1 ? new MathMap<>() {
					@Override
					public NFE evaluate(Real t) {
						return getInteger(t.round());
					}
				} : null;
				this.realEmbeddings.add(new EmbeddedNumberField<>(this, r, alpha.realPart(), round));
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
				this.complexEmbeddings.add(new EmbeddedNumberField<>(this, c, alpha, round));
			}
		}
		embeddingsComputed = true;
	}

	public List<EmbeddedNumberField<Real>> realEmbeddings() {
		computeEmbeddings();
		return Collections.unmodifiableList(this.realEmbeddings);
	}

	public List<EmbeddedNumberField<ComplexNumber>> complexEmbeddings() {
		computeEmbeddings();
		return Collections.unmodifiableList(this.complexEmbeddings);
	}

	public List<EmbeddedNumberField<Ext<PAdicNumber>>> padicEmbeddings(PAdicField base) {
		IntE prime = base.getPrime();
		UnivariatePolynomialRing<PAdicNumber> padicPolynomials = base.getUnivariatePolynomialRing();
		Integers z = Integers.z();
		DiscreteValuationRing<Fraction, PFE> localized = z.localize(prime);
		PrimeField fp = PrimeField.getPrimeField(prime);
		Extension<PFE, PFE, FFE, FiniteField> trivialExtension = fp
				.getExtension(fp.getUnivariatePolynomialRing().getVar());
		TheMontesResult<Fraction, PFE, PFE, FFE, FiniteField> theMontes = localized
				.theMontesAlgorithm(minimalPolynomial(), trivialExtension);
		List<EmbeddedNumberField<Ext<PAdicNumber>>> result = new ArrayList<>();
		for (OkutsuType<Fraction, PFE, PFE, FFE, FiniteField> type : theMontes.getTypes()) {
			type = localized.singleFactorLifting(type, base.getAccuracy());
			CompleteDVRExtension<PAdicNumber, PFE, PFE, FFE, FiniteField> embeddingField = CompleteDVRExtension
					.getCompleteDVRExtension(
							padicPolynomials.getEmbedding(type.representative(), base.fromRationalMap()), base,
							trivialExtension);
			Ext<PAdicNumber> alpha = embeddingField.alpha();
			result.add(new EmbeddedNumberField<>(this, embeddingField, alpha, new MathMap<>() {
				@Override
				public NFE evaluate(Ext<PAdicNumber> t) {
					UnivariatePolynomial<Fraction> asPolynomial = q.getUnivariatePolynomialRing()
							.getEmbedding(t.asPolynomial(), base.toRationalMap());
					UnivariatePolynomial<NFE> asNFEPolynomial = getUnivariatePolynomialRing().getEmbedding(asPolynomial,
							getEmbeddingMap());
					return getUnivariatePolynomialRing().evaluate(asNFEPolynomial, alpha());
				}
			}));
		}
		return result;
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
		for (EmbeddedNumberField<Real> realEmbedding : realEmbeddings()) {
			embedding.add(realEmbedding.embedding(t));
		}
		for (EmbeddedNumberField<ComplexNumber> complexEmbedding : complexEmbeddings()) {
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
		for (EmbeddedNumberField<Real> realEmbedding : realEmbeddings()) {
			result.add(r.log(realEmbedding.value(t)));
		}
		for (EmbeddedNumberField<ComplexNumber> complexEmbedding : complexEmbeddings()) {
			result.add(r.multiply(2, r.log(complexEmbedding.value(t))));
		}
		return new Vector<>(result);
	}

	@Override
	public FactorizationResult<Polynomial<NFE>, NFE> factorization(UnivariatePolynomial<NFE> t) {
		UnivariatePolynomialRing<NFE> ring = getUnivariatePolynomialRing();
		SortedMap<Polynomial<NFE>, Integer> result = new TreeMap<>();
		NFE unit = t.leadingCoefficient();
		t = removeDenominators(t);
		FactorizationResult<Polynomial<NFE>, NFE> squareFree = ring.squareFreeFactorization(t);
		for (Polynomial<NFE> squareFreeFactor : squareFree.primeFactors()) {
			for (Polynomial<NFE> factor : factorizeSquareFree(removeDenominators(squareFreeFactor))) {
				result.put(factor, squareFree.multiplicity(squareFreeFactor));
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
			Vector<Fraction> c = asVector(p.univariateCoefficient(i));
			for (int j = 0; j < degree(); j++) {
				lcm = z.lcm(lcm, c.get(j + 1).getDenominator());
			}
		}
		return ring.multiply(getInteger(lcm), p);
	}

	private List<Polynomial<NFE>> factorizeSquareFree(UnivariatePolynomial<NFE> t) {
		UnivariatePolynomialRing<NFE> ring = getUnivariatePolynomialRing();
		if (t.degree() == 1) {
			return Collections.singletonList(ring.normalize(t));
		}
		Rationals q = Rationals.q();
		Integers z = Integers.z();
		UnivariatePolynomialRing<Fraction> rationalRing = q.getUnivariatePolynomialRing();
		GenericAlgebraicRingExtension<NFE> cr = new GenericAlgebraicRingExtension<>(t, this);
		Iterator<UnivariatePolynomial<IntE>> it = z.getUnivariatePolynomialRing().polynomials(degree() - 1);
		while (true) {
			NFE probeTerm = fromSmallDegreePolynomial(rationalRing.getEmbedding(it.next(), new MathMap<>() {
				@Override
				public Fraction evaluate(IntE t) {
					return q.getEmbedding(t);
				}
			}));
			GenericAlgebraicExtensionElement<NFE> probe = cr.add(cr.alpha(),
					cr.getEmbedding(multiply(probeTerm, alpha())));
			List<Vector<Fraction>> generators = new ArrayList<>();
			GenericAlgebraicExtensionElement<NFE> power = cr.one();
			for (int i = 0; i < t.degree() * degree(); i++) {
				generators.add(asRationalVector(power, t.degree()));
				power = cr.multiply(power, probe);
			}
			FiniteVectorSpace<Fraction> asVectorSpace = new FiniteVectorSpace<>(q, t.degree() * degree());
			Matrix<Fraction> m = Matrix.fromColumns(generators);
			MatrixAlgebra<Fraction> algebra = asVectorSpace.matrixAlgebra();
			if (algebra.rank(m) != t.degree() * degree()) {
				continue;
			}
			UnivariatePolynomial<Fraction> x = rationalRing.getPolynomial(
					algebra.solve(m, asRationalVector(cr.fromPolynomial(ring.getVar()), t.degree())).asList());
			UnivariatePolynomial<Fraction> gamma = rationalRing
					.getPolynomial(algebra.solve(m, asRationalVector(cr.getEmbedding(alpha()), t.degree())).asList());
			Vector<Fraction> lastPower = asRationalVector(power, t.degree());
			UnivariatePolynomial<Fraction> minimalPolynomial = rationalRing
					.getPolynomial(algebra.solve(m, lastPower).asList());
			minimalPolynomial = rationalRing.toUnivariate(
					rationalRing.subtract(rationalRing.getVarPower(t.degree() * degree()), minimalPolynomial));
			FactorizationResult<Polynomial<Fraction>, Fraction> rationalFactors = q.factorization(minimalPolynomial);
			if (rationalFactors.primeFactors().size() == 1) {
				return Collections.singletonList(ring.normalize(t));
			}
			List<Polynomial<NFE>> factors = new ArrayList<>();
			for (Polynomial<Fraction> factor : rationalFactors.primeFactors()) {
				NumberField extension = new NumberField(rationalRing.toUnivariate(factor));
				if (!extension.minimalPolynomial().equals(factor)) {
					throw new ArithmeticException("Something messed with the minimalpolynomial!");
				}
				UnivariatePolynomialRing<NFE> extensionRing = extension.getUnivariatePolynomialRing();
				NFE embeddedX = extensionRing.evaluate(extensionRing.getEmbedding(x, extension.getEmbeddingMap()),
						extension.alpha());
				NFE embeddedGamma = extensionRing
						.evaluate(extensionRing.getEmbedding(gamma, extension.getEmbeddingMap()), extension.alpha());
				List<Vector<Fraction>> gammaXBase = new ArrayList<>();
				NFE xPower = extension.one();
				for (int i = 0; i < extension.degree() / degree(); i++) {
					NFE gammaPower = extension.one();
					for (int j = 0; j < degree(); j++) {
						gammaXBase.add(extension.asVector(extension.multiply(xPower, gammaPower)));
						gammaPower = extension.multiply(gammaPower, embeddedGamma);
					}
					xPower = extension.multiply(xPower, embeddedX);
				}
				Matrix<Fraction> gammaXMatrix = Matrix.fromColumns(gammaXBase);
				UnivariatePolynomial<NFE> nfeFactor = minimalPolynomial(embeddedX, extension.degree() / degree(),
						extension, this, new MathMap<>() {

							@Override
							public Vector<NFE> evaluate(NFE t) {
								Vector<Fraction> overGammaX = extension.matrixAlgebra().solve(gammaXMatrix,
										extension.asVector(t));
								List<NFE> asList = new ArrayList<>();
								for (int i = 0; i < extension.degree() / degree(); i++) {
									asList.add(fromPolynomial(rationalRing.getPolynomial(
											overGammaX.asList().subList(i * degree(), (i + 1) * degree()))));
								}
								return new Vector<>(asList);
							}
						});
				factors.add(ring.normalize(nfeFactor));
			}
			return factors;
		}
	}

	private Vector<Fraction> asRationalVector(GenericAlgebraicExtensionElement<NFE> e, int dimension) {
		UnivariatePolynomialRing<NFE> ring = getUnivariatePolynomialRing();
		List<NFE> asNFEList = ring.asVector(e.asPolynomial(), dimension - 1).asList();
		List<Fraction> asList = new ArrayList<>();
		for (NFE nfe : asNFEList) {
			asList.addAll(asVector(nfe).asList());
		}
		return new Vector<>(asList);
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
		IdealClassGroup icg = idealClassGroup();
		UnivariatePolynomialRing<NFE> polynomialRing = getUnivariatePolynomialRing();
		UnivariatePolynomial<NFE> unity = polynomialRing.toUnivariate(polynomialRing
				.subtract(polynomialRing.getVarPower(classNumber().intValueExact()), polynomialRing.one()));
		FactorizationResult<Polynomial<NFE>, NFE> rootsOfUnity = factorization(unity);
		for (Polynomial<NFE> factor : rootsOfUnity.primeFactors()) {
			if (factor.degree() != 1) {
				throw new UnsupportedOperationException("Need more roots of unity!");
			}
		}
		List<NFE> unitGroupGenerators = maximalOrder().unitGroupGenerators();
		List<NFE> generators = new ArrayList<>();
		generators.addAll(unitGroupGenerators);
		List<NFE> alpha = new ArrayList<>();
		List<IdealClass> diagonalGenerators = icg.getDiagonalModuleGenerators();
		List<IntE> ranks = icg.getDiagonalRanks();
		for (int i = 0; i < diagonalGenerators.size(); i++) {
			IdealClass ic = diagonalGenerators.get(i);
			NumberFieldIdeal power = maximalOrder().power(ic.representative(), ranks.get(i).intValueExact());
			if (!power.isPrincipal()) {
				throw new ArithmeticException("Ideal Class Group wrong!");
			}
			alpha.add(power.generators().get(0));
		}
		generators.addAll(alpha);
		return null;
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
		// TODO Auto-generated method stub
		return null;
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
