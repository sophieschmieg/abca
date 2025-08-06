package varieties.curves.elliptic;

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;

import fields.exceptions.InfinityException;
import fields.helper.AbstractAlgebra;
import fields.helper.AbstractElement;
import fields.integers.FiniteRationalVectorSpace;
import fields.integers.Integers;
import fields.integers.Integers.IntE;
import fields.integers.Integers.IntegerIdeal;
import fields.integers.Rationals;
import fields.integers.Rationals.Fraction;
import fields.interfaces.AlgebraicExtensionElement;
import fields.interfaces.Element;
import fields.interfaces.Field;
import fields.interfaces.FieldExtension;
import fields.interfaces.Ideal;
import fields.interfaces.UnivariatePolynomialRing;
import fields.interfaces.Field.Extension;
import fields.numberfields.NumberField;
import fields.numberfields.NumberField.NFE;
import fields.numberfields.NumberFieldOrder;
import fields.quaternions.RationalQuaternionOrder;
import fields.quaternions.RationalQuaternions;
import fields.vectors.FreeModule;
import fields.vectors.Matrix;
import fields.vectors.MatrixAlgebra;
import fields.vectors.Vector;
import util.MiscAlgorithms;
import varieties.curves.elliptic.EllipticCurve.EllipticCurveExtensionResult;
import varieties.curves.elliptic.EllipticCurve.Isomorphism;
import varieties.curves.elliptic.EndomorphismRing.Endomorphism;
import varieties.projective.ProjectiveMorphism;
import varieties.projective.ProjectivePoint;

public class EndomorphismRing<T extends Element<T>> extends AbstractAlgebra<IntE, Endomorphism<T>> {

	public static class Endomorphism<T extends Element<T>> extends AbstractElement<Endomorphism<T>> {
		private Vector<IntE> asVector;
		private Isogeny<T> asIsogeny;
		private EndomorphismRing<T> endomorphismRing;

		private Endomorphism(Vector<IntE> asVector, EndomorphismRing<T> endomorphismRing) {
			this.asVector = asVector;
			this.endomorphismRing = endomorphismRing;
			this.asIsogeny = new Isogeny<>() {

				@Override
				public int compareTo(Isogeny<T> o) {
					throw new ArithmeticException("not comparable");
				}

				@Override
				public BigInteger getDegree() {
					return endomorphismRing.reducedNorm(Endomorphism.this).getValue();
				}

				@Override
				public EllipticCurve<T> getDomain() {
					return endomorphismRing.curve;
				}

				@Override
				public EllipticCurve<T> getRange() {
					return endomorphismRing.curve;
				}

				@Override
				public ProjectivePoint<T> evaluate(ProjectivePoint<T> point) {
					ProjectivePoint<T> result = endomorphismRing.curve.neutral();
					for (int i = 0; i < asVector.dimension(); i++) {
						result = endomorphismRing.curve.add(endomorphismRing.curve.multiply(asVector.get(i + 1),
								endomorphismRing.basis.get(i).evaluate(point)), result);
					}
					return result;
				}

				@Override
				public <B extends Element<B>, E extends AlgebraicExtensionElement<B, E>, FE extends FieldExtension<B, E, FE>> ProjectivePoint<E> evaluateOverExtension(
						ProjectivePoint<E> point, Extension<T, B, E, FE> extension) {
					EllipticCurve<E> extendedCurve = endomorphismRing.curve.extendBaseField(extension).getCurve();
					ProjectivePoint<E> result = extendedCurve.neutral();
					for (int i = 0; i < asVector.dimension(); i++) {
						result = extendedCurve.add(extendedCurve.multiply(asVector.get(i + 1),
								endomorphismRing.basis.get(i).evaluateOverExtension(point, extension)), result);
					}
					return result;
				}

				@Override
				public ProjectiveMorphism<T> asMorphism() {
					// TODO Auto-generated method stub
					return null;
				}

				@Override
				public Isogeny<T> getDual() {
					return endomorphismRing.getDual(Endomorphism.this).asIsogeny;
				}

				@Override
				public List<ProjectivePoint<T>> kernelGenerators() {
					// TODO Auto-generated method stub
					return null;
				}

			};
		}

		@Override
		public int compareTo(Endomorphism<T> o) {
			return asVector.compareTo(o.asVector);
		}

		public ProjectivePoint<T> evaluate(ProjectivePoint<T> t) {
			return asIsogeny.evaluate(t);
		}
	}

	private EllipticCurve<T> curve;
	private List<Isogeny<T>> basis;
	private Integers z;
	private NumberFieldOrder asNumberFieldOrder;
	private RationalQuaternionOrder asQuaternionOrder;
	private FreeModule<IntE> asFreeModule;

	EndomorphismRing(EllipticCurve<T> curve) {
		this.curve = curve;
		this.z = Integers.z();
		Rationals q = Rationals.q();
		Field<T> field = curve.getField();
		Isogeny<T> identity = curve.identity();
		basis = new ArrayList<>();
		if (!curve.getField().isFinite()) {
			throw new UnsupportedOperationException("Only curves over finite fields currently supported");
		}

		Isogeny<T> frobenius = new Frobenius<>(curve, curve.minimalDefinedExtension());
		if (!frobenius.getRange().equals(curve)) {
			frobenius = new CompositionIsogeny<>(frobenius,
					frobenius.getRange().getIsomorphisms(curve).iterator().next());
		}
		if (curve.isSupersingular()) {
			if (curve.jInvariant().equals(curve.getField().getInteger(1728))) {
				RationalQuaternions quaternions = RationalQuaternions.quaternions(q.getInteger(-1),
						q.getInteger(curve.getField().characteristic().negate()));
				this.asQuaternionOrder = new RationalQuaternionOrder(quaternions, quaternions.i(),
						quaternions.divide(quaternions.add(quaternions.i(), quaternions.j()),
								quaternions.getInteger(2)),
						quaternions.divide(quaternions.add(quaternions.one(), quaternions.k()),
								quaternions.getInteger(2)));
				if (!field.hasSqrt(field.getInteger(-1))) {
					throw new ArithmeticException("Extend field for endomorphism ring!");
				}
				T sqrt = field.sqrt(field.getInteger(-1)).keySet().iterator().next();
				Isomorphism<T> isomorphism = new Isomorphism<>(curve, sqrt, field.zero(), field.zero(), field.zero());
				Isogeny<T> traceZero = new AdditionIsogeny<>(isomorphism, frobenius);
				Isogeny<T> traceOne = new AdditionIsogeny<>(curve.identity(),
						new CompositionIsogeny<>(frobenius, isomorphism));
				if (!isDivisible(traceZero, BigInteger.TWO)) {
					traceZero = new AdditionIsogeny<>(isomorphism, new CompositionIsogeny<>(isomorphism, frobenius));
					traceOne = new AdditionIsogeny<>(curve.identity(), frobenius);
				}
				basis.add(identity);
				basis.add(isomorphism);
				basis.add(DivisionIsogeny.getDivisionIsogeny(traceZero, field.characteristic().add(BigInteger.ONE), BigInteger.TWO));
				basis.add(DivisionIsogeny.getDivisionIsogeny(traceOne, field.characteristic().add(BigInteger.ONE), BigInteger.TWO));
			} else {
			}
		} else {
			if (curve.jInvariant().equals(field.getInteger(1728))) {
				asNumberFieldOrder = NumberField
						.getNumberField(q.getUnivariatePolynomialRing().getPolynomial(q.one(), q.zero(), q.one()))
						.maximalOrder().asOrder();
				T sqrt = field.sqrt(field.getInteger(-1)).keySet().iterator().next();
				Isomorphism<T> isomorphism = new Isomorphism<>(curve, sqrt, field.zero(), field.zero(), field.zero());
				basis.add(identity);
				basis.add(isomorphism);
			} else if (curve.jInvariant().equals(field.zero())) {
				asNumberFieldOrder = NumberField
						.getNumberField(q.getUnivariatePolynomialRing().getPolynomial(q.one(), q.one(), q.one()))
						.maximalOrder().asOrder();
				T cbrt = field
						.roots(field.getInteger(asNumberFieldOrder.numberField()
								.trace(asNumberFieldOrder.getModuleGenerators().get(1)).asInteger()), 3)
						.keySet().iterator().next();
				Isomorphism<T> isomorphism = new Isomorphism<>(curve, cbrt, field.zero(), field.zero(), field.zero());
				basis.add(identity);
				basis.add(isomorphism);
			} else {
				IntE trace = z.getInteger(curve.trace());
				IntE norm = z.getInteger(frobenius.getDegree());
				IntE disc = z.subtract(z.multiply(4, norm), z.multiply(trace, trace));
				FactorizationResult<IntE, IntE> discFactors = z.uniqueFactorization(disc);
				IntE discriminant = z.one();
				Fraction real = q.getFraction(trace, z.getInteger(2));
				Fraction imaginary = q.getFraction(1, 2);
				IntE conductor = z.one();
				for (IntE prime : discFactors.primeFactors()) {
					discriminant = z.multiply(z.power(prime, discFactors.multiplicity(prime) % 2), discriminant);
					imaginary = q.multiply(
							z.power(prime, MiscAlgorithms.DivRoundDown(discFactors.multiplicity(prime), 2)), imaginary);
					conductor = z.multiply(
							z.power(prime,
									MiscAlgorithms.DivRoundDown(discFactors.multiplicity(prime), 2) - height(prime)),
							conductor);
				}
				if (!z.remainder(discriminant, z.getInteger(4)).equals(z.getInteger(3))) {
					conductor = z.divideChecked(conductor, z.getInteger(2));
				}
				UnivariatePolynomialRing<Fraction> rationalPolynomials = q.getUnivariatePolynomialRing();
				NumberField nf = NumberField.getNumberField(
						rationalPolynomials.getPolynomial(q.getInteger(discriminant), q.zero(), q.one()));
				asNumberFieldOrder = nf
						.getOrder(nf.multiply(conductor, nf.maximalOrder().getModuleGenerators().get(1)));
				NFE frobeniousNFE = nf.add(nf.getFraction(real), nf.scalarMultiply(imaginary, nf.alpha()));
				List<Vector<Fraction>> inFrobeniousBasis = new ArrayList<>();
				inFrobeniousBasis.add(Vector.mapVector(q.getEmbeddingMap(), asNumberFieldOrder.asVector(nf.one())));
				inFrobeniousBasis
						.add(Vector.mapVector(q.getEmbeddingMap(), asNumberFieldOrder.asVector(frobeniousNFE)));
				MatrixAlgebra<Fraction> matrixAlgebra = new FiniteRationalVectorSpace(2).matrixAlgebra();
				Matrix<Fraction> baseChange = matrixAlgebra.inverse(Matrix.fromColumns(inFrobeniousBasis));
				Vector<Fraction> nonTrivial = baseChange.column(2);
				IntE lcm = z.lcm(nonTrivial.get(1).getDenominator(), nonTrivial.get(2).getDenominator());
				basis.add(curve.identity());
				basis.add(DivisionIsogeny.getDivisionIsogeny(
						new AdditionIsogeny<>(
								new MultiplicationIsogeny<>(q.multiply(lcm, nonTrivial.get(1)).asInteger().getValue(),
										curve),
								new CompositionIsogeny<>(frobenius,
										new MultiplicationIsogeny<>(
												q.multiply(lcm, nonTrivial.get(2)).asInteger().getValue(), curve))),
						nf.norm(asNumberFieldOrder.getModuleGenerators().get(1)).asInteger().getValue()
								.multiply(lcm.getValue()),
						lcm.getValue()));
			}
		}
		for (Isogeny<T> isogeny : basis) {
			System.out.println(isogeny + " norm: " + reducedNorm(isogeny) + " trace: " + reducedTrace(isogeny));
		}
		asFreeModule = new FreeModule<>(z, basis.size());
	}

	private int height(IntE prime) {
		List<ProjectivePoint<T>> basis = curve.getTorsionPointBasis(prime);
		if (basis.size() < 2) {
			return 0;
		}
		List<Isogeny<T>> isogenies = new ArrayList<>();
		isogenies.add(new KernelPointIsogeny<>(curve, basis.get(0), prime.intValueExact()));
		isogenies.add(new KernelPointIsogeny<>(curve, basis.get(1), prime.intValueExact()));
		isogenies.add(new KernelPointIsogeny<>(curve, curve.add(basis.get(0), basis.get(1)), prime.intValueExact()));
		int height = 0;
		while (true) {
			height++;
			List<Isogeny<T>> newIsogenies = new ArrayList<>();
			for (Isogeny<T> isogeny : isogenies) {
				List<ProjectivePoint<T>> torsionBasis = isogeny.getRange().getTorsionPointBasis(prime);
				if (torsionBasis.size() < 2) {
					return height;
				}
				Isogeny<T> next = null;
				do {
					ProjectivePoint<T> kernelPoint = isogeny.getRange().getRandomTorsionPoint(prime);
					if (kernelPoint.equals(isogeny.getRange().neutral())) {
						next = null;
						continue;
					}
					next = new KernelPointIsogeny<>(isogeny.getRange(), kernelPoint, prime.intValueExact());
				} while (next == null || next.getRange().equals(isogeny.getDomain()));
				newIsogenies.add(next);
			}
			isogenies = newIsogenies;
		}
	}

	private boolean isDivisible(Isogeny<T> isogeny, BigInteger divisor) {
		List<ProjectivePoint<T>> divisorBasis = curve.getTorsionPointBasis(divisor);
		if (divisorBasis.size() != 2) {
			throw new ArithmeticException("Not enough divisors to divide!");
		}
		for (ProjectivePoint<T> point : divisorBasis) {
			if (!isogeny.evaluate(point).equals(curve.neutral())) {
				return false;
			}
		}
		return true;
	}

	@Override
	public Endomorphism<T> getEmbedding(IntE t) {
		List<IntE> list = new ArrayList<>();
		list.add(t);
		for (int i = 1; i < basis.size(); i++) {
			list.add(z.zero());
		}
		return fromVector(new Vector<>(list));
	}

	public IntE reducedNorm(Isogeny<T> isogeny) {
		return z.getInteger(isogeny.getDegree());
	}

	public IntE reducedTrace(Isogeny<T> isogeny) {
		return z.getInteger(curve.asMultiplicationIsogenyModOrder(new AdditionIsogeny<>(isogeny, isogeny.getDual())));
	}

	public IntE reducedNorm(Endomorphism<T> t) {
		if (asQuaternionOrder != null) {
			return asQuaternionOrder.getQuaternions().reducedNorm(asQuaternionOrder.fromVector(t.asVector)).asInteger();
		}
		return asNumberFieldOrder.numberField().norm(asNumberFieldOrder.fromVector(t.asVector)).asInteger();
	}

	public IntE reducedTrace(Endomorphism<T> t) {
		if (asQuaternionOrder != null) {
			return asQuaternionOrder.getQuaternions().reducedTrace(asQuaternionOrder.fromVector(t.asVector))
					.asInteger();
		}
		return asNumberFieldOrder.numberField().trace(asNumberFieldOrder.fromVector(t.asVector)).asInteger();
	}

	@Override
	public Endomorphism<T> fromVector(Vector<IntE> vector) {
		return new Endomorphism<>(vector, this);
	}

	public Endomorphism<T> getDual(Endomorphism<T> t) {
		if (asQuaternionOrder != null) {
			return fromVector(asQuaternionOrder
					.asVector(asQuaternionOrder.getQuaternions().conjugate(asQuaternionOrder.fromVector(t.asVector))));
		}
		NFE nfe = asNumberFieldOrder.fromVector(t.asVector);
		List<NFE> conjugates = asNumberFieldOrder.numberField().conjugates(nfe);
		for (NFE conjugate : conjugates) {
			if (!conjugate.equals(nfe)) {
				return fromVector(asNumberFieldOrder.asVector(conjugate));
			}
		}
		throw new ArithmeticException("No conjugates!");
	}

	@Override
	public boolean isGeneratingAlgebra(List<Endomorphism<T>> s) {
		// TODO Auto-generated method stub
		return false;
	}

	@Override
	public List<Endomorphism<T>> getAlgebraGenerators() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public Integers getRing() {
		return Integers.z();
	}

	@Override
	public Endomorphism<T> zero() {
		return fromVector(asFreeModule.zero());
	}

	@Override
	public Endomorphism<T> add(Endomorphism<T> s1, Endomorphism<T> s2) {
		return fromVector(asFreeModule.add(s1.asVector, s2.asVector));
	}

	@Override
	public Endomorphism<T> negative(Endomorphism<T> s) {
		return fromVector(asFreeModule.negative(s.asVector));
	}

	@Override
	public boolean isFree() {
		return true;
	}

	@Override
	public IntegerIdeal annihilator() {
		return z.getZeroIdeal();
	}

	@Override
	public List<Vector<IntE>> getSyzygies() {
		return Collections.emptyList();
	}

	@Override
	public boolean isLinearIndependent(List<Endomorphism<T>> s) {
		// TODO Auto-generated method stub
		return false;
	}

	@Override
	public boolean isGeneratingModule(List<Endomorphism<T>> s) {
		// TODO Auto-generated method stub
		return false;
	}

	@Override
	public List<Vector<IntE>> nonTrivialCombinations(List<Endomorphism<T>> s) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public List<Endomorphism<T>> getModuleGenerators() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public Vector<IntE> asVector(Endomorphism<T> s) {
		return s.asVector;
	}

	@Override
	public Exactness exactness() {
		return Exactness.EXACT;
	}

	@Override
	public Endomorphism<T> getRandomElement() {
		return fromVector(asFreeModule.getRandomElement());
	}

	@Override
	public boolean isFinite() {
		return false;
	}

	@Override
	public BigInteger getNumberOfElements() throws InfinityException {
		return asFreeModule.getNumberOfElements();
	}

	@Override
	public Iterator<Endomorphism<T>> iterator() {
		throw new InfinityException();
	}

	@Override
	public Endomorphism<T> one() {
		return fromVector(asFreeModule.getUnitVector(1));
	}

	@Override
	public BigInteger characteristic() {
		return BigInteger.ZERO;
	}

	@Override
	public Endomorphism<T> multiply(Endomorphism<T> t1, Endomorphism<T> t2) {
		if (asQuaternionOrder != null) {
			return fromVector(asQuaternionOrder.asVector(asQuaternionOrder
					.multiply(asQuaternionOrder.fromVector(t1.asVector), asQuaternionOrder.fromVector(t2.asVector))));
		}
		return fromVector(asNumberFieldOrder.asVector(asNumberFieldOrder
				.multiply(asNumberFieldOrder.fromVector(t1.asVector), asNumberFieldOrder.fromVector(t2.asVector))));
	}

	@Override
	public boolean isUnit(Endomorphism<T> t) {
		if (asQuaternionOrder != null) {
			return asQuaternionOrder.isUnit(asQuaternionOrder.fromVector(t.asVector));
		}
		return asNumberFieldOrder.isUnit(asNumberFieldOrder.fromVector(t.asVector));
	}

	@Override
	public Endomorphism<T> inverse(Endomorphism<T> t) {
		if (asQuaternionOrder != null) {
			return fromVector(
					asQuaternionOrder.asVector(asQuaternionOrder.inverse(asQuaternionOrder.fromVector(t.asVector))));
		}
		return fromVector(
				asNumberFieldOrder.asVector(asNumberFieldOrder.inverse(asNumberFieldOrder.fromVector(t.asVector))));
	}

	@Override
	public boolean isCommutative() {
		return asQuaternionOrder == null;
	}

	@Override
	public boolean isIntegral() {
		return true;
	}

	@Override
	public boolean isReduced() {
		return true;
	}

	@Override
	public boolean isIrreducible() {
		return true;
	}

	@Override
	public boolean isZeroDivisor(Endomorphism<T> t) {
		return t.equals(zero());
	}

	@Override
	public boolean isEuclidean() {
		return false;
	}

	@Override
	public boolean isUniqueFactorizationDomain() {
		return false;
	}

	@Override
	public FactorizationResult<Endomorphism<T>, Endomorphism<T>> uniqueFactorization(Endomorphism<T> t) {
		throw new ArithmeticException("not a UFD");
	}

	@Override
	public boolean isPrincipalIdealDomain() {
		return false;
	}

	@Override
	public boolean isDedekindDomain() {
		return false;
	}

	@Override
	public boolean isDivisible(Endomorphism<T> dividend, Endomorphism<T> divisor) {
		// TODO Auto-generated method stub
		return false;
	}

	@Override
	public QuotientAndRemainderResult<Endomorphism<T>> quotientAndRemainder(Endomorphism<T> dividend,
			Endomorphism<T> divisor) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public BigInteger euclidMeasure(Endomorphism<T> t) {
		return null;
	}

	@Override
	public Endomorphism<T> projectToUnit(Endomorphism<T> t) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public Iterable<Endomorphism<T>> getUnits() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public int krullDimension() {
		// TODO Auto-generated method stub
		return 0;
	}

	@Override
	public IdealResult<Endomorphism<T>, ?> getIdealWithTransforms(List<Endomorphism<T>> generators) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public Ideal<Endomorphism<T>> intersect(Ideal<Endomorphism<T>> t1, Ideal<Endomorphism<T>> t2) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public Ideal<Endomorphism<T>> radical(Ideal<Endomorphism<T>> t) {
		// TODO Auto-generated method stub
		return null;
	}
}
