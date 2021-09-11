package fields.numberfields;

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.SortedMap;
import java.util.TreeMap;

import fields.floatingpoint.Complex;
import fields.floatingpoint.Complex.ComplexNumber;
import fields.floatingpoint.Reals;
import fields.floatingpoint.Reals.Real;
import fields.helper.AbstractElement;
import fields.helper.AbstractFieldExtension;
import fields.helper.AbstractInnerProductSpace;
import fields.helper.GaloisGroup;
import fields.helper.GenericAlgebraicRingExtension;
import fields.helper.GenericAlgebraicRingExtension.GenericAlgebraicExtensionElement;
import fields.integers.Integers;
import fields.integers.Integers.IntE;
import fields.integers.Rationals;
import fields.integers.Rationals.Fraction;
import fields.interfaces.AlgebraicExtensionElement;
import fields.interfaces.InnerProductSpace;
import fields.interfaces.MathMap;
import fields.interfaces.Polynomial;
import fields.interfaces.UnivariatePolynomial;
import fields.interfaces.UnivariatePolynomialRing;
import fields.interfaces.ValueField;
import fields.numberfields.NumberField.NFE;
import fields.vectors.FiniteVectorSpace;
import fields.vectors.Matrix;
import fields.vectors.MatrixAlgebra;
import fields.vectors.Vector;
import util.MiscAlgorithms;

public class NumberField extends AbstractFieldExtension<Fraction, NFE, NumberField>
		implements InnerProductSpace<Fraction, NFE> {
	private Rationals q;
	private Reals r;
	private Complex c;
	private boolean embeddingsComputed;
	private List<EmbeddedNumberField<Real>> realEmbeddings;
	private List<EmbeddedNumberField<ComplexNumber>> complexEmbeddings;
	private NumberFieldIntegers maximalOrder;
	private IdealGroup idealGroup;
	private IdealClassGroup idealClassGroup;
	private List<Fraction> basisNorms;
	private static final int PRECISION_FOR_EMBEDDINGS = 1024;

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
	}

	public NumberFieldIntegers maximalOrder() {
		if (maximalOrder == null) {
			maximalOrder = new NumberFieldIntegers(this);
		}
		return maximalOrder;
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
			if (alpha.complexPart().equals(r.zero()) || alpha.complexPart().exponent() < -r.precision()) {
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
//
//	public List<EmbeddedNumberField<Ext<PAdicNumber>>> padicEmbeddings(BigInteger prime, int accuracy) {
//		NumberFieldIntegers r = maximalOrder();
//		PAdicField base = new PAdicField(prime, accuracy);
//		Integers z = Integers.z();
//		z.reduceUnivariatePolynomial(minimalPolynomial(), prime);
//		UnivariatePolynomial<PAdicNumber> minPoly = base.getUnivariatePolynomialRing().getEmbedding(minimalPolynomial(),
//				new MathMap<>() {
//					@Override
//					public PAdicNumber evaluate(Fraction t) {
//						return base.fromRational(t);
//					}
//				});
//		FactorizationResult<Polynomial<PFE>> factors = base.reduction()
//				.factorization(base.ringOfIntegers().reduceUnivariatePolynomial(minPoly));
//		List<EmbeddedNumberField<Ext<PAdicNumber, PFE, FFE>>> result = new ArrayList<>();
//		for (Polynomial<PFE> factor : factors.getFactors().keySet()) {
//			LocalExtension<PAdicNumber, PAdicNumber, PFE, FFE> zp = base.getExtension(base.ringOfIntegers()
//					.henselLiftFactor(minPoly, base.reduction().getUnivariatePolynomialRing().toUnivariate(factor)));
//			Ext<PAdicNumber, PFE, FFE> gamma = zp
//					.roots(zp.getUnivariatePolynomialRing().getEmbedding(minPoly, zp.getPrimeEmbeddingMap())).keySet()
//					.iterator().next();
//			result.add(new EmbeddedNumberField<>(this, zp, gamma, new MathMap<>() {
//				@Override
//				public NFE evaluate(Ext<PAdicNumber, PFE, FFE> t) {
//					Polynomial<PAdicNumber> asPoly = zp.asGenericPrimeExtensionFieldElement(t).asPolynomial();
//					Polynomial<Fraction> roundedPoly = Rationals.q().getUnivariatePolynomialRing().getEmbedding(asPoly,
//							new MathMap<>() {
//								@Override
//								public Fraction evaluate(PAdicNumber t) {
//									return Rationals.q().getEmbedding(base.roundToInteger(t, accuracy));
//								}
//							});
//					return getUnivariatePolynomialRing().evaluate(
//							getUnivariatePolynomialRing().getEmbedding(roundedPoly, getPrimeEmbeddingMap()), gamma());
//				}
//			}));
//		}
//		return result;
//	}

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

	public BigInteger classNumber() {
		return idealClassGroup().getNumberOfElements();
//		Integers z = Integers.z();
//		if (degree() > 2) {
//			throw new UnsupportedOperationException("Class number not implemented for degree > 2!");
//		}
//		if (degree() == 1) {
//			return BigInteger.ONE;
//		}
//		IntE d = discriminant();
//		if (d.compareTo(z.zero()) > 0) {
//			throw new UnsupportedOperationException("Class number not implemented for real quadratic fields!");
//		}
//		FactorizationResult<IntE> primeFactors = z.uniqueFactorization(d);
//		primeLoop: for (IntE prime : primeFactors.primeFactors()) {
//			for (int i = 0; i < z.divide(d, prime).intValueExact(); i++) {
//				IntE a = z.add(z.multiply(i, prime), z.one());
//				if (!z.gcd(d, a).equals(z.one())) {
//					continue;
//				}
//				if (MiscAlgorithms.kroneckerSymbol(d.getValue(), a.getValue()) != 1) {
//					continue primeLoop;
//				}
//			}
//			System.out.println("Prime: " + prime);
//			System.out.println("Mod 4: " + prime.getValue().mod(BigInteger.valueOf(4)));
//		}
//		return null;
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

	@Override
	public ValueField<Fraction> getValueField() {
		return Rationals.q().withInfValue();
	}

	@Override
	public Fraction innerProduct(NFE s1, NFE s2) {
		if (basisNorms == null) {
			basisNorms = new ArrayList<>();
			for (int i = 0; i < degree(); i++) {
				basisNorms.add(q.withInfValue().abs(q.power(minimalPolynomial().univariateCoefficient(0), i)));
			}
		}
		Vector<Fraction> t1 = asVector(s1);
		Vector<Fraction> t2 = asVector(s2);
		Fraction result = q.zero();
		for (int i = 0; i < degree(); i++) {
			result = q.add(result, q.multiply(basisNorms.get(i), t1.get(i + 1), t2.get(i + 1)));
		}
		return result;
	}

	@Override
	public List<NFE> latticeReduction(List<NFE> s) {
		return AbstractInnerProductSpace.latticeReduction(this, new MathMap<>() {
			@Override
			public Fraction evaluate(Fraction t) {
				return q.withInfValue().round(t);
			}
		}, s);
	}

}
