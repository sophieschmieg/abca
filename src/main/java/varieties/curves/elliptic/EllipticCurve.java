package varieties.curves.elliptic;

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Deque;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;

import fields.finitefields.FiniteField;
import fields.finitefields.FiniteField.FFE;
import fields.finitefields.PrimeField;
import fields.finitefields.PrimeField.PFE;
import fields.floatingpoint.Complex;
import fields.floatingpoint.Complex.ComplexNumber;
import fields.floatingpoint.Reals;
import fields.helper.AbstractElement;
import fields.helper.FieldEmbedding;
import fields.integers.Integers;
import fields.integers.Integers.IntE;
import fields.integers.Rationals;
import fields.integers.Rationals.Fraction;
import fields.interfaces.DiscreteValuationField;
import fields.interfaces.Element;
import fields.interfaces.Field;
import fields.interfaces.Group;
import fields.interfaces.Ideal;
import fields.interfaces.MathMap;
import fields.interfaces.Polynomial;
import fields.interfaces.PolynomialRing;
import fields.interfaces.Ring.ExtendedEuclideanResult;
import fields.interfaces.Ring.FactorizationResult;
import fields.interfaces.Ring.QuotientAndRemainderResult;
import fields.interfaces.UnivariatePolynomial;
import fields.interfaces.UnivariatePolynomialRing;
import fields.local.Value;
import fields.numberfields.EmbeddedNumberField;
import fields.numberfields.FractionalIdeal;
import fields.numberfields.LocalizedNumberField;
import fields.numberfields.NumberField;
import fields.numberfields.NumberField.NFE;
import fields.numberfields.NumberFieldIntegers.NumberFieldIdeal;
import fields.numberfields.NumberFieldOrder;
import fields.numberfields.PicardGroup.OrderIdealClass;
import fields.polynomials.AbstractPolynomialRing;
import fields.polynomials.CoordinateRing;
import fields.polynomials.CoordinateRing.CoordinateRingElement;
import fields.polynomials.Monomial;
import fields.polynomials.PolynomialIdeal;
import fields.vectors.FreeSubModule;
import fields.vectors.Matrix;
import fields.vectors.Vector;
import util.MiscAlgorithms;
import util.Pair;
import util.graph.Graph;
import util.graph.Graph.GraphBuilder;
import varieties.FunctionField;
import varieties.RationalFunction;
import varieties.Scheme;
import varieties.affine.AffinePoint;
import varieties.curves.SmoothCurve;
import varieties.curves.WeilDivisor;
import varieties.projective.AbstractProjectiveScheme;
import varieties.projective.GenericProjectiveScheme;
import varieties.projective.ProjectiveMorphism;
import varieties.projective.ProjectivePoint;

public class EllipticCurve<T extends Element<T>> extends AbstractProjectiveScheme<T>
		implements SmoothCurve<T>, Group<ProjectivePoint<T>>, Element<EllipticCurve<T>> {
//	public static class Isomorphism<T extends Element<T>> extends AbstractElement<Isogeny<T>> implements Isogeny<T> {
//		private final EllipticCurve<T> domain;
//		private final EllipticCurve<T> range;
//		private final Field<T> field;
//		private final T u;
//
//		private Isomorphism(EllipticCurve<T> domain, T u) {
//			this.u = u;
//			this.domain = domain;
//			this.field = domain.getField();
//			this.range = new EllipticCurve<T>(field, field.multiply(field.power(u, 4), domain.getA()),
//					field.multiply(field.power(u, 6), domain.getB()));
//			range.setNumberOfPointsFrom(domain, this);
//		}
//
//		@Override
//		public String toString() {
//			return "[" + u + "]";
//		}
//
//		@Override
//		public int compareTo(Isogeny<T> o) {
//			if (o instanceof Isomorphism<?>) {
//				Isomorphism<T> other = (Isomorphism<T>) o;
//				return u.compareTo(other.u);
//			}
//			return getClass().getName().compareTo(o.getClass().getName());
//		}
//
//		@Override
//		public EllipticCurve<T> getDomain() {
//			return domain;
//		}
//
//		@Override
//		public EllipticCurve<T> getRange() {
//			return range;
//		}
//
//		@Override
//		public ProjectivePoint<T> evaluate(ProjectivePoint<T> point) {
//			T x = point.getCoord(1);
//			T y = point.getCoord(2);
//			T z = point.getCoord(3);
//			return new ProjectivePoint<>(field, field.multiply(field.power(u, 2), x),
//					field.multiply(field.power(u, 3), y), z);
//		}
//
//		@Override
//		public ProjectiveMorphism<T> asMorphism() {
//			PolynomialRing<T> projectiveRing = domain.asGenericProjectiveScheme().homogenousPolynomialRing();
//			List<Polynomial<T>> asPolynomials = new ArrayList<>();
//			asPolynomials.add(projectiveRing.multiply(field.power(u, 2), projectiveRing.getVar(1)));
//			asPolynomials.add(projectiveRing.multiply(field.power(u, 3), projectiveRing.getVar(2)));
//			asPolynomials.add(projectiveRing.getVar(3));
//			return new ProjectiveMorphism<>(domain.asGenericProjectiveScheme(), range.asGenericProjectiveScheme(),
//					asPolynomials);
//		}
//
//		@Override
//		public BigInteger getDegree() {
//			return BigInteger.ONE;
//		}
//
//		@Override
//		public Isogeny<T> getDual() {
//			return new Isomorphism<>(range, field.inverse(u));
//		}
//
////		@Override
//		public List<ProjectivePoint<T>> kernelGenerators() {
//			return Collections.emptyList();
//		}
//	}

	public static class IdentityIsomorphism<T extends Element<T>> extends AbstractElement<Isogeny<T>>
			implements Isogeny<T> {
		private final EllipticCurve<T> domain;

		private IdentityIsomorphism(EllipticCurve<T> domain) {
			this.domain = domain;
		}

		@Override
		public String toString() {
			return "id";
		}

		@Override
		public int compareTo(Isogeny<T> o) {
			if (o instanceof IdentityIsomorphism<?>) {
				IdentityIsomorphism<T> other = (IdentityIsomorphism<T>) o;
				return domain.compareTo(other.domain);
			}
			return getClass().getName().compareTo(o.getClass().getName());
		}

		@Override
		public EllipticCurve<T> getDomain() {
			return domain;
		}

		@Override
		public EllipticCurve<T> getRange() {
			return domain;
		}

		@Override
		public ProjectivePoint<T> evaluate(ProjectivePoint<T> point) {
			return point;
		}

		@Override
		public ProjectiveMorphism<T> asMorphism() {
			PolynomialRing<T> projectiveRing = domain.asGenericProjectiveScheme().homogenousPolynomialRing();
			List<Polynomial<T>> asPolynomials = new ArrayList<>();
			asPolynomials.add(projectiveRing.getVar(1));
			asPolynomials.add(projectiveRing.getVar(2));
			asPolynomials.add(projectiveRing.getVar(3));
			return new ProjectiveMorphism<>(domain.asGenericProjectiveScheme(), domain.asGenericProjectiveScheme(),
					asPolynomials);
		}

		@Override
		public BigInteger getDegree() {
			return BigInteger.ONE;
		}

		@Override
		public Isogeny<T> getDual() {
			return this;
		}

//		@Override
		public List<ProjectivePoint<T>> kernelGenerators() {
			return Collections.emptyList();
		}
	}

	public static class Isomorphism<T extends Element<T>> extends AbstractElement<Isogeny<T>> implements Isogeny<T> {
		private final EllipticCurve<T> domain;
		private final EllipticCurve<T> range;
		private final Field<T> field;
		private final T scale;
		private final T xOffset;
		private final T xyOffset;
		private final T yOffset;

		private Isomorphism(EllipticCurve<T> domain, T scale, T xOffset, T xyOffset, T yOffset) {
			this.domain = domain;
			this.field = domain.getField();
			this.scale = scale;
			this.xOffset = xOffset;
			this.xyOffset = xyOffset;
			this.yOffset = yOffset;
			PolynomialRing<T> r = domain.affineRing;
			List<Polynomial<T>> eval = new ArrayList<>();
			eval.add(r.add(r.getEmbedding(field.power(scale, 2), new int[] { 1, 0 }), r.getEmbedding(xOffset)));
			eval.add(r.add(r.getEmbedding(field.power(scale, 3), new int[] { 0, 1 }),
					r.getEmbedding(xyOffset, new int[] { 1, 0 }), r.getEmbedding(yOffset)));
			Polynomial<T> substituted = r.substitute(domain.definingPolynomial, eval);
			this.range = fromScaledPolynomial(substituted);
		}

		@Override
		public String toString() {
			Isomorphism<T> dual = getDual();
			return "(" + field.power(dual.scale, 2) + "*X + " + dual.xOffset + ", " + field.power(dual.scale, 3)
					+ "*Y + " + dual.xyOffset + "X + " + dual.yOffset + ")";
		}

		@Override
		public int compareTo(Isogeny<T> o) {
			if (o instanceof Isomorphism<?>) {
				Isomorphism<T> other = (Isomorphism<T>) o;
				int cmp = scale.compareTo(other.scale);
				if (cmp != 0) {
					return cmp;
				}
				cmp = xOffset.compareTo(other.xOffset);
				if (cmp != 0) {
					return cmp;
				}
				cmp = xyOffset.compareTo(other.xyOffset);
				if (cmp != 0) {
					return cmp;
				}
				cmp = yOffset.compareTo(other.yOffset);
				if (cmp != 0) {
					return cmp;
				}
				return domain.compareTo(other.domain);
			}
			return getClass().getName().compareTo(o.getClass().getName());
		}

		@Override
		public EllipticCurve<T> getDomain() {
			return domain;
		}

		@Override
		public EllipticCurve<T> getRange() {
			return range;
		}

		@Override
		public ProjectivePoint<T> evaluate(ProjectivePoint<T> point) {
			if (point.getCoord(3).equals(field.zero())) {
				return point;
			}
			T x = point.getDehomogenisedCoord(1, 3);
			T y = point.getDehomogenisedCoord(2, 3);
			T xPrime = field.divide(field.subtract(x, xOffset), field.power(scale, 2));
			T yPrime = field.divide(field.subtract(y, field.add(field.multiply(xyOffset, xPrime), yOffset)),
					field.power(scale, 3));
			return new ProjectivePoint<>(field, xPrime, yPrime, field.one());
		}

		@Override
		public ProjectiveMorphism<T> asMorphism() {
			Isomorphism<T> dual = getDual();
			PolynomialRing<T> projectiveRing = domain.asGenericProjectiveScheme().homogenousPolynomialRing();
			List<Polynomial<T>> asPolynomials = new ArrayList<>();
			asPolynomials.add(
					projectiveRing.add(projectiveRing.multiply(field.power(dual.scale, 2), projectiveRing.getVar(1)),
							projectiveRing.getEmbedding(dual.xOffset, new int[] { 0, 0, 1 })));
			asPolynomials.add(
					projectiveRing.add(projectiveRing.multiply(field.power(dual.scale, 3), projectiveRing.getVar(2)),
							projectiveRing.getEmbedding(dual.xyOffset, new int[] { 1, 0, 0 }),
							projectiveRing.getEmbedding(dual.yOffset, new int[] { 0, 0, 1 })));
			asPolynomials.add(projectiveRing.getVar(3));
			return new ProjectiveMorphism<>(domain.asGenericProjectiveScheme(), range.asGenericProjectiveScheme(),
					asPolynomials);
		}

		@Override
		public BigInteger getDegree() {
			return BigInteger.ONE;
		}

		@Override
		public Isomorphism<T> getDual() {
			T inverseScale = field.inverse(scale);
			T inverseXOffset = field.multiply(-1, xOffset, field.power(inverseScale, 2));
			T inverseXYOffset = field.multiply(-1, xyOffset, field.power(inverseScale, 5));
			T inverseYOffset = field.add(
					field.multiply(field.multiply(xOffset, xyOffset), field.power(inverseScale, 5)),
					field.multiply(-1, yOffset, field.power(inverseScale, 3)));
			return new Isomorphism<>(range, inverseScale, inverseXOffset, inverseXYOffset, inverseYOffset);
		}

//		@Override
		public List<ProjectivePoint<T>> kernelGenerators() {
			return Collections.emptyList();
		}
	}

	private Field<T> field;
	private T a1;
	private T a2;
	private T a3;
	private T a4;
	private T a6;
	private T b2;
	private T b4;
	private T b6;
	private T b8;
	private T c4;
	private T c6;
	private T discriminant;
	private T j;
	private ProjectivePoint<T> pointAtInfinity;
	private PolynomialRing<T> projectiveRing;
	private PolynomialRing<T> affineRing;
	private UnivariatePolynomialRing<T> univariateRing;
	private CoordinateRing<T> coordinateRing;
	private Polynomial<T> rhs;
	private Polynomial<T> lhs;
	private Map<Integer, Polynomial<T>> divisionPolynomials;
	private Polynomial<T> definingPolynomial;
	private BigInteger numberOfPoints = BigInteger.valueOf(-1);
	private boolean superSingular;
	private boolean superSingularDeterminied;
	private Map<IntE, List<ProjectivePoint<T>>> torsionPointBasis;
	private Isogeny<T> weierstrassForm;

	private static <T extends Element<T>> GenericProjectiveScheme<T> asGenericProjectiveScheme(Field<T> field, T a1,
			T a2, T a3, T a4, T a6) {
		PolynomialRing<T> ring = AbstractPolynomialRing.getPolynomialRing(field, 3, Monomial.GREVLEX);
		Polynomial<T> f = ring.getEmbedding(field.one(), new int[] { 3, 0, 0 });
		f = ring.add(f, ring.getEmbedding(field.negative(field.one()), new int[] { 0, 2, 1 }));
		f = ring.add(f, ring.getEmbedding(field.negative(a1), new int[] { 1, 1, 1 }));
		f = ring.add(f, ring.getEmbedding(field.negative(a3), new int[] { 0, 1, 2 }));
		f = ring.add(f, ring.getEmbedding(a2, new int[] { 2, 0, 1 }));
		f = ring.add(f, ring.getEmbedding(a4, new int[] { 1, 0, 2 }));
		f = ring.add(f, ring.getEmbedding(a6, new int[] { 0, 0, 3 }));
		return new GenericProjectiveScheme<>(field, ring, Collections.singletonList(f));
	}

	public static <T extends Element<T>> EllipticCurve<T> fromScaledPolynomial(Polynomial<T> polynomial) {
		if (polynomial.degree() != 3) {
			throw new ArithmeticException("degree unequal to three!");
		}
		PolynomialRing<T> r = polynomial.getPolynomialRing();
		T zero = r.getRing().zero();
		T one = r.getRing().one();
		T negativeOne = r.getRing().negative(one);
		if (r.numberOfVariables() != 2) {
			throw new ArithmeticException("more than two variables!");
		}
		if (!polynomial.coefficient(r.getMonomial(new int[] { 0, 3 })).equals(zero)
				|| !polynomial.coefficient(r.getMonomial(new int[] { 1, 2 })).equals(zero)
				|| !polynomial.coefficient(r.getMonomial(new int[] { 2, 1 })).equals(zero)) {
			throw new ArithmeticException("Y^3, XY^2, X^2Y have non zero coefficient!");
		}
		T scale = polynomial.coefficient(r.getMonomial(new int[] { 3, 0 }));
		if (scale.equals(zero)) {
			throw new ArithmeticException("X^3 has coefficient zero!");
		}
		polynomial = r.divideScalar(polynomial, scale);
		if (!polynomial.coefficient(r.getMonomial(new int[] { 0, 2 })).equals(negativeOne)) {
			throw new ArithmeticException("Y^2 have coefficient unequal to the negative of the coefficient of X^3!");
		}
		return new EllipticCurve<>((Field<T>) r.getRing(),
				r.getRing().negative(polynomial.coefficient(r.getMonomial(new int[] { 1, 1 }))),
				polynomial.coefficient(r.getMonomial(new int[] { 2, 0 })),
				r.getRing().negative(polynomial.coefficient(r.getMonomial(new int[] { 0, 1 }))),
				polynomial.coefficient(r.getMonomial(new int[] { 1, 0 })),
				polynomial.coefficient(r.getMonomial(new int[] { 0, 0 })));
	}

	public static class FromPolynomialResult<T extends Element<T>> {
		private List<Polynomial<T>> substituions;
		private EllipticCurve<T> curve;
		private MathMap<ProjectivePoint<T>, ProjectivePoint<T>> map;

		private FromPolynomialResult(EllipticCurve<T> curve, List<Polynomial<T>> substituions) {
			this.curve = curve;
			this.substituions = substituions;
			this.map = new MathMap<>() {

				@Override
				public ProjectivePoint<T> evaluate(ProjectivePoint<T> t) {
					if (t.getCoord(3).equals(curve.field.zero())) {
						return curve.pointAtInfinity;
					}
					List<T> points = new ArrayList<>();
					for (Polynomial<T> sub : substituions) {
						points.add(sub.getPolynomialRing().evaluate(sub, t.getDehomogenous(3)));
					}
					points.add(curve.field.one());
					return new ProjectivePoint<>(curve.field, points);
				}
			};
		}

		public List<Polynomial<T>> getSubstituions() {
			return substituions;
		}

		public EllipticCurve<T> getCurve() {
			return curve;
		}

		public MathMap<ProjectivePoint<T>, ProjectivePoint<T>> getMap() {
			return map;
		}
	}

	public static <T extends Element<T>> FromPolynomialResult<T> fromPolynomial(Field<T> field,
			Polynomial<T> polynomial) {
		PolynomialRing<T> ring = AbstractPolynomialRing.getPolynomialRing(field, 2, Monomial.REVLEX);
		if (polynomial.numberOfVariables() == 3 && polynomial.getPolynomialRing().isHomogeneous(polynomial)) {
			polynomial = ring.getEmbedding(polynomial.getPolynomialRing().dehomogenize(polynomial, 3),
					new int[] { 0, 1 });
		} else if (polynomial.numberOfVariables() == 2) {
			polynomial = ring.getEmbedding(polynomial, new int[] { 0, 1 });
		} else {
			throw new ArithmeticException("Unexpected number of variables!");
		}
		if (polynomial.degree(1) != 3 || polynomial.degree(2) != 2) {
			throw new ArithmeticException("Degree wrong!");
		}
		T scale = field.divide(polynomial.coefficient(ring.getMonomial(new int[] { 0, 2 })),
				polynomial.coefficient(ring.getMonomial(new int[] { 3, 0 })));
		List<Polynomial<T>> substitutions = new ArrayList<>();
		substitutions.add(ring.divideScalar(ring.getVar(1), scale));
		substitutions.add(ring.divideScalar(ring.getVar(2), scale));
		List<Polynomial<T>> scaling = new ArrayList<>();
		scaling.add(ring.multiply(scale, ring.getVar(1)));
		scaling.add(ring.multiply(scale, ring.getVar(2)));
		EllipticCurve<T> curve = fromScaledPolynomial(ring.substitute(polynomial, scaling));
		return new FromPolynomialResult<>(curve, substitutions);
	}

	private static <T extends Element<T>> Pair<T, T> fromJInvariantCoefficients(Field<T> field, T j) {
		if (j.equals(field.zero())) {
			return new Pair<>(field.zero(), field.one());
		}
		if (j.equals(field.getInteger(1728))) {
			return new Pair<>(field.negative(field.one()), field.zero());
		}
		T a = field.divide(field.multiply(4, field.subtract(field.getInteger(1728), j)), field.multiply(27, j));
		T b = field.power(a, 2);
		return new Pair<>(a, b);
	}

	public static <T extends Element<T>> EllipticCurve<T> fromJInvariant(Field<T> field, T j) {
		if (field.characteristic().equals(BigInteger.TWO)) {
			if (j.equals(field.zero())) {
				return new EllipticCurve<>(field, field.zero(), field.zero(), field.one(), field.zero(), field.zero());
			}
			return new EllipticCurve<>(field, field.one(), field.zero(), field.zero(), field.zero(),
					field.negative(field.inverse(j)));
		} else if (field.characteristic().equals(BigInteger.valueOf(3))) {
			if (j.equals(field.zero())) {
				return new EllipticCurve<>(field, field.zero(), field.zero(), field.zero(), field.negative(field.one()),
						field.zero());
			}

		}
		Pair<T, T> coefficients = fromJInvariantCoefficients(field, j);
		return new EllipticCurve<>(field, coefficients.getFirst(), coefficients.getSecond());
	}

	public static <T extends Element<T>> EllipticCurve<T> fromLengendre(Field<T> field, T lambda) {
		PolynomialRing<T> polynomialRing = AbstractPolynomialRing.getPolynomialRing(field, 2, Monomial.REVLEX);
		Polynomial<T> lhs = polynomialRing.getVarPower(2, 2);
		Polynomial<T> rhs = polynomialRing.multiply(polynomialRing.getVar(1),
				polynomialRing.subtract(polynomialRing.getVar(1), polynomialRing.one()),
				polynomialRing.subtract(polynomialRing.getVar(1), polynomialRing.getEmbedding(lambda)));
		return fromScaledPolynomial(polynomialRing.subtract(rhs, lhs));
	}

	public static class FromTauResult {
		private ComplexNumber tau;
		private EllipticCurve<ComplexNumber> curve;
		private MathMap<ComplexNumber, ProjectivePoint<ComplexNumber>> weierstrassP;
		private MathMap<ComplexNumber, ProjectivePoint<ComplexNumber>> weierstrassPDerivative;

		public FromTauResult(ComplexNumber tau, EllipticCurve<ComplexNumber> curve,
				MathMap<ComplexNumber, ProjectivePoint<ComplexNumber>> weierstrassP,
				MathMap<ComplexNumber, ProjectivePoint<ComplexNumber>> weierstrassPDerivative) {
			super();
			this.tau = tau;
			this.curve = curve;
			this.weierstrassP = weierstrassP;
			this.weierstrassPDerivative = weierstrassPDerivative;
		}

		public ComplexNumber getTau() {
			return tau;
		}

		public EllipticCurve<ComplexNumber> getCurve() {
			return curve;
		}

		public MathMap<ComplexNumber, ProjectivePoint<ComplexNumber>> getWeierstrassP() {
			return weierstrassP;
		}

		public MathMap<ComplexNumber, ProjectivePoint<ComplexNumber>> getWeierstrassPDerivative() {
			return weierstrassPDerivative;
		}
	}

	private static IntE eisensteinSeries(IntE n, int weight) {
		Integers z = Integers.z();
		IntE result = z.zero();
		for (IntE factor : z.factors(n)) {
			result = z.add(z.power(factor, weight), result);
		}
		return result;
	}

	private static List<IntE> zagierCoefficients = new ArrayList<>();

	private static IntE zagierCoefficient(int n) {
		Integers z = Integers.z();
		if (n < -1 || n % 4 == 1 || n % 4 == 2) {
			return z.zero();
		} else if (n == -1) {
			return z.getInteger(-1);
		} else if (n == 0) {
			return z.getInteger(2);
		}
		int index = 2 * MiscAlgorithms.DivRoundUp(n, 4) - 1 - (n % 2);
		if (zagierCoefficients.size() <= index) {
			int div4 = MiscAlgorithms.DivRoundUp(n, 4);
			IntE zagier3 = z.multiply(-480, eisensteinSeries(z.getInteger(div4), 3));
			for (int r = 2; r * r - 1 <= 4 * div4; r++) {
				zagier3 = z.subtract(zagier3, z.multiply(2 * r * r, zagierCoefficient(4 * div4 - r * r)));
			}
			IntE zagier4 = z.negative(zagier3);
			zagier3 = z.divideChecked(zagier3, z.getInteger(2));
			for (int r = 2; r * r - 1 <= 4 * div4; r++) {
				zagier4 = z.subtract(zagier4, z.multiply(2, zagierCoefficient(4 * div4 - r * r)));
			}
			if (zagierCoefficients.size() != 2 * div4 - 2) {
				throw new ArithmeticException("Zagier recursion wrong!");
			}
			zagierCoefficients.add(zagier3);
			zagierCoefficients.add(zagier4);
		}
		return zagierCoefficients.get(index);
	}

	private static IntE jInvariantCoefficient(int n) {
		Integers z = Integers.z();
		if (n == 0) {
			return z.getInteger(744);
		}
		IntE result = z.zero();
		if (n % 4 != 0) {
			result = zagierCoefficient(n);
			for (int r = 1; r * r - 1 <= n; r++) {
				result = z.add(z.multiply(2, zagierCoefficient(n - r * r)), result);
			}
		}
		for (int i = 0; 4 * i * i + 4 * i <= 16 * n; i++) {
			int r = 2 * i + 1;
			result = z.add(z.subtract(z.multiply(n % 2 == 0 ? 1 : -1, zagierCoefficient(4 * n - r * r)),
					zagierCoefficient(16 * n - r * r)), result);
		}
		return z.divideChecked(result, z.getInteger(n));
	}

	private static ComplexNumber jInvariant(Complex c, ComplexNumber tau) {
		c = c.withPrecision(2 * c.precision());
		Reals r = c.getReals();
		if (tau.complexPart().compareTo(r.zero()) < 0) {
			tau = c.negative(tau);
		}
		tau = c.subtract(tau, c.getInteger(tau.realPart().round()));
		while (c.value(tau).compareTo(r.one()) < 0) {
			tau = c.negative(c.inverse(tau));
			tau = c.subtract(tau, c.getInteger(tau.realPart().round()));
		}
		ComplexNumber q = c.exp(c.multiply(2, c.pi(), c.i(), tau));
		ComplexNumber j = c.zero();
		ComplexNumber prev;
		int i = -1;
		do {
			prev = j;
			j = c.add(c.multiply(jInvariantCoefficient(i), c.power(q, i)), j);
			i++;
		} while (!c.close(j, prev));
		return c.withPrecision(c.precision() / 2).getEmbedding(j);
	}

	private static ComplexNumber jInvariant(Complex c, NumberField field, OrderIdealClass ideal) {
		return jInvariant(c, field, ideal.representative().getNumerator().asSubModule());
	}

	private static ComplexNumber jInvariant(Complex c, NumberField field, FreeSubModule<IntE, NFE> submodule) {
		EmbeddedNumberField<ComplexNumber, Complex> embedding = field.complexEmbeddings(c).get(0);
		List<NFE> generators = submodule.getModuleGenerators();
		ComplexNumber tau = c.divide(embedding.embedding(generators.get(0)), embedding.embedding(generators.get(1)));
		return jInvariant(c, tau);
	}

	private static MathMap<ComplexNumber, ProjectivePoint<ComplexNumber>> weierstrassP(Complex c, ComplexNumber tau) {
		return new MathMap<>() {

			@Override
			public ProjectivePoint<ComplexNumber> evaluate(ComplexNumber t) {
				if (t.equals(c.zero())) {
					return new ProjectivePoint<>(c, c.one(), c.zero());
				}
				ComplexNumber result = c.power(c.inverse(t), 2);
				ComplexNumber prev;
				int i = 1;
				do {
					prev = result;
					for (int n = 0; n <= i; n++) {
						int m = i - n;
						for (int j = 0; j < 4; j++) {
							ComplexNumber lambda;
							switch (j) {
							case 0:
								lambda = c.add(c.getInteger(m), c.multiply(n, tau));
								break;
							case 1:
								lambda = c.add(c.getInteger(m), c.multiply(-n, tau));
								break;
							case 2:
								lambda = c.add(c.getInteger(-m), c.multiply(n, tau));
								break;
							case 3:
								lambda = c.add(c.getInteger(-m), c.multiply(-n, tau));
								break;
							default:
								lambda = null;
								break;
							}
							ComplexNumber summand = c.power(c.inverse(c.subtract(t, lambda)), 2);
							summand = c.subtract(summand, c.power(c.inverse(lambda), 2));
							result = c.add(result, summand);
						}
					}
				} while (!c.close(prev, result));
				return new ProjectivePoint<>(c, result, c.one());
			}
		};
	}

	private static MathMap<ComplexNumber, ProjectivePoint<ComplexNumber>> weierstrassPDerivative(Complex c,
			ComplexNumber tau) {
		return new MathMap<>() {

			@Override
			public ProjectivePoint<ComplexNumber> evaluate(ComplexNumber t) {
				if (t.equals(c.zero())) {
					return new ProjectivePoint<>(c, c.one(), c.zero());
				}
				ComplexNumber result = c.power(c.inverse(t), 3);
				ComplexNumber prev;
				int i = 1;
				do {
					prev = result;
					for (int n = 0; n <= i; n++) {
						int m = i - n;
						for (int j = 0; j < 4; j++) {
							ComplexNumber lambda;
							switch (j) {
							case 0:
								lambda = c.add(c.getInteger(m), c.multiply(n, tau));
								break;
							case 1:
								lambda = c.add(c.getInteger(m), c.multiply(-n, tau));
								break;
							case 2:
								lambda = c.add(c.getInteger(-m), c.multiply(n, tau));
								break;
							case 3:
								lambda = c.add(c.getInteger(-m), c.multiply(-n, tau));
								break;
							default:
								lambda = null;
								break;
							}
							ComplexNumber summand = c.power(c.inverse(c.subtract(t, lambda)), 3);
							result = c.add(result, summand);
						}
					}
				} while (!c.close(prev, result));
				return new ProjectivePoint<>(c, c.multiply(-2, result), c.one());
			}
		};
	}

	public static FromTauResult fromTau(Complex c, ComplexNumber tau) {
		ComplexNumber j = jInvariant(c, tau);
		return new FromTauResult(tau, fromJInvariant(c, j), weierstrassP(c, tau), weierstrassPDerivative(c, tau));
	}

	public static class WithComplexMultiplication {
		private EllipticCurve<NFE> curve;
		private Isogeny<NFE> endomorphism;
		private NFE tau;
		private NumberField imaginaryQuadraticExtension;
		private FieldEmbedding<Fraction, NFE, NumberField> rayClassField;

		private WithComplexMultiplication(EllipticCurve<NFE> curve, Isogeny<NFE> endomorphism, NFE tau,
				NumberField imaginaryQuadraticExtension, FieldEmbedding<Fraction, NFE, NumberField> rayClassField) {
			this.curve = curve;
			this.endomorphism = endomorphism;
			this.tau = tau;
			this.imaginaryQuadraticExtension = imaginaryQuadraticExtension;
			this.rayClassField = rayClassField;
		}

		public EllipticCurve<NFE> getCurve() {
			return curve;
		}

		public Isogeny<NFE> getEndomorphism() {
			return endomorphism;
		}

		public NFE getTau() {
			return tau;
		}

		public NumberField getImaginaryQuadraticExtension() {
			return imaginaryQuadraticExtension;
		}

		public FieldEmbedding<Fraction, NFE, NumberField> getRayClassField() {
			return rayClassField;
		}
	}

	public static class ClassFieldResult {
		private FieldEmbedding<Fraction, NFE, NumberField> classField;
		private NFE jInvariant;
		private UnivariatePolynomial<NFE> minimalPolynomial;
		private boolean rationalMinimalPolynomial;
		private UnivariatePolynomial<NFE> jInvariantAsPolynomial;

		private ClassFieldResult(FieldEmbedding<Fraction, NFE, NumberField> classField, NFE jInvariant,
				UnivariatePolynomial<NFE> minimalPolynomial, boolean rationalMinimalPolynomial,
				UnivariatePolynomial<NFE> jInvariantAsPolynomial) {
			this.classField = classField;
			this.jInvariant = jInvariant;
			this.minimalPolynomial = minimalPolynomial;
			this.rationalMinimalPolynomial = rationalMinimalPolynomial;
			this.jInvariantAsPolynomial = jInvariantAsPolynomial;
		}

		public FieldEmbedding<Fraction, NFE, NumberField> getClassField() {
			return classField;
		}

		public NFE getjInvariant() {
			return jInvariant;
		}

		public UnivariatePolynomial<NFE> getMinimalPolynomial() {
			return minimalPolynomial;
		}

		public boolean isRationalMinimalPolynomial() {
			return rationalMinimalPolynomial;
		}

		public UnivariatePolynomial<NFE> getjInvariantAsPolynomial() {
			return jInvariantAsPolynomial;
		}
	}

	public static ClassFieldResult computeHilbertClassField(NumberField field) {
		return computeRayClassField(field.maximalOrder().asOrder());
	}

	public static ClassFieldResult computeRayClassField(NumberFieldOrder order) {
		NumberField field = order.numberField();
		Integers z = Integers.z();
		if (field.degree() != 2 || field.discriminant().compareTo(z.zero()) > 0) {
			throw new ArithmeticException("Not an imaginary quadratic number field");
		}
		Complex c = Complex.c(MiscAlgorithms.roundUpToPowerOfTwo(256 * order.classNumber().intValueExact()));
		List<ComplexNumber> conjugates = new ArrayList<>();
		for (OrderIdealClass ic : order.picardGroup()) {
			conjugates.add(jInvariant(c, field, ic));
		}
		Pair<UnivariatePolynomial<NFE>, UnivariatePolynomial<NFE>> simplifiedMinimalPolynomial = field
				.simplifyMinimalPolynomial(field.integerMinimalPolynomial(c, conjugates));
		UnivariatePolynomial<NFE> minimalPolynomial = simplifiedMinimalPolynomial.getFirst();
		FieldEmbedding<Fraction, NFE, NumberField> rayClassField = field.getEmbeddedExtension(minimalPolynomial);
		NFE jInvariant = rayClassField.fromPolynomial(simplifiedMinimalPolynomial.getSecond());
		boolean rational = true;
		for (int i = 0; i <= minimalPolynomial.degree(); i++) {
			if (minimalPolynomial.univariateCoefficient(i).asPolynomial().degree() > 0) {
				rational = false;
				break;
			}
		}
		if (rational) {
			return new ClassFieldResult(rayClassField, jInvariant, minimalPolynomial, true,
					simplifiedMinimalPolynomial.getSecond());
		}
		NumberField rcf = rayClassField.getField();
		for (NFE t : rcf.maximalOrder().getModuleGenerators()) {
			UnivariatePolynomial<Fraction> mipo = rcf.minimalPolynomial(t);
			if (rcf.degree() == mipo.degree()) {
				NumberField simplerRCF = NumberField.getNumberField(mipo);
				NFE power = rcf.one();
				List<Vector<Fraction>> asVectors = new ArrayList<>();
				for (int i = 0; i < rcf.degree(); i++) {
					asVectors.add(rcf.asVector(power));
					power = rcf.multiply(power, t);
				}
				Matrix<Fraction> baseChange = Matrix.fromColumns(asVectors);
				Vector<Fraction> alpha = rcf.matrixAlgebra().solve(baseChange,
						rcf.asVector(rayClassField.getEmbeddedAlpha()));
				rayClassField = new FieldEmbedding<>(simplerRCF, field, simplerRCF.fromVector(alpha));
				jInvariant = simplerRCF.fromVector(rcf.matrixAlgebra().solve(baseChange, rcf.asVector(jInvariant)));
				break;
			}
		}
		return new ClassFieldResult(rayClassField, jInvariant, rayClassField.minimalPolynomial(), false,
				rayClassField.asPolynomial(jInvariant));
	}

	public static WithComplexMultiplication withComplexMultiplication(NumberField field, NFE tau) {
		Integers z = Integers.z();
		Rationals q = Rationals.q();
		if (field.degree() != 2 || field.discriminant().compareTo(z.zero()) > 0) {
			throw new ArithmeticException("Not an imaginary quadratic number field");
		}
		Fraction minusOne = q.negative(q.one());
		Fraction one = q.one();
		Fraction trace = field.trace(tau);
		Fraction norm = field.norm(tau);
		while (trace.compareTo(minusOne) < 0 || trace.compareTo(one) > 0 || norm.compareTo(q.one()) < 0) {
			trace = field.trace(tau);
			if (trace.compareTo(minusOne) < 0 || trace.compareTo(one) > 0) {
				IntE rounded = q.divide(trace, q.getInteger(2)).round();
				tau = field.subtract(tau, field.getEmbedding(rounded));
			}
			norm = field.norm(tau);
			if (norm.compareTo(q.one()) < 0) {
				tau = field.negative(field.inverse(tau));
			}
		}
		NFE orderGenerator = field.multiply(field.maximalOrder().getNumerator(tau),
				field.maximalOrder().getDenominator(tau));
		NumberFieldOrder order = field.getOrder(orderGenerator);
		ClassFieldResult classFieldResult = computeRayClassField(order);
		FieldEmbedding<Fraction, NFE, NumberField> rayClassField = classFieldResult.getClassField();
		NFE jInvariant = classFieldResult.getjInvariant();
		UnivariatePolynomial<NFE> minimalPolynomial = classFieldResult.getMinimalPolynomial();
		boolean rational = classFieldResult.isRationalMinimalPolynomial();
		NumberField classField = rayClassField.getField();
		if (rational) {
			MathMap<NFE, Fraction> retraction = new MathMap<>() {
				@Override
				public Fraction evaluate(NFE t) {
					return t.asPolynomial().univariateCoefficient(0);
				}
			};
			UnivariatePolynomial<Fraction> rationalMinimalPolynomial = q.getUnivariatePolynomialRing()
					.getEmbedding(minimalPolynomial, retraction);
			classField = NumberField.getNumberField(rationalMinimalPolynomial);
			jInvariant = classField.fromPolynomial(q.getUnivariatePolynomialRing()
					.getEmbedding(classFieldResult.getjInvariantAsPolynomial(), retraction));
		}
		Pair<NFE, NFE> coefficients = fromJInvariantCoefficients(classField, jInvariant);
		NFE a = coefficients.getFirst();
		NFE b = coefficients.getSecond();
		Set<Ideal<NFE>> factors = new TreeSet<>();
		if (!a.equals(classField.zero())) {
			FractionalIdeal aIdeal = FractionalIdeal.principalIdeal(classField.maximalOrder(), a);
			factors.addAll(classField.maximalOrder().idealFactorization(aIdeal.getNumerator()).primeFactors());
			factors.addAll(classField.maximalOrder().idealFactorization(aIdeal.getDenominator()).primeFactors());
		}
		if (!b.equals(classField.zero())) {
			FractionalIdeal bIdeal = FractionalIdeal.principalIdeal(classField.maximalOrder(), b);
			factors.addAll(classField.maximalOrder().idealFactorization(bIdeal.getNumerator()).primeFactors());
			factors.addAll(classField.maximalOrder().idealFactorization(bIdeal.getDenominator()).primeFactors());
		}
		NFE multiplier = classField.one();
		for (Ideal<NFE> ideal : factors) {
			NumberFieldIdeal prime = (NumberFieldIdeal) ideal;
			LocalizedNumberField localized = classField.maximalOrder().localizeAndQuotient(ideal);
			NFE uniformizer = localized.uniformizer();
			if (prime.isPrincipal()) {
				uniformizer = prime.principalGenerator();
			} else if (prime.type().ramificationIndex() > 1) {
				uniformizer = classField.getInteger(prime.intGenerator());
			}
			int uniformizerValue = localized.valuation(uniformizer).value();
			Value aValue = localized.valuation(a);
			Value bValue = localized.valuation(b);
			Value exponent = Value.INFINITY;
			if (!aValue.isInfinite()) {
				exponent = exponent.min(new Value(MiscAlgorithms.DivRoundDown(aValue.value(), 2 * uniformizerValue)));
			}
			if (!bValue.isInfinite()) {
				exponent = exponent.min(new Value(MiscAlgorithms.DivRoundDown(bValue.value(), 3 * uniformizerValue)));
			}
			multiplier = classField.multiply(classField.power(uniformizer, -exponent.value()), multiplier);
		}
		a = classField.multiply(a, multiplier, multiplier);
		b = classField.multiply(b, classField.power(multiplier, 3));
		System.err.println(jInvariant);
		EllipticCurve<NFE> curve = new EllipticCurve<>(classField, a, b);
		System.err.println(curve.jInvariant());
		return new WithComplexMultiplication(curve, null, tau, field, rayClassField);
	}

	public static enum KodairaSymbol {
		I0, I1, I2, Iv, mIv, II, III, IV, I0star, Ivstar, IVstar, IIIstar, IIstar;
	}

	private static class CurveParameters<T extends Element<T>> {
		private T a1;
		private T a2;
		private T a3;
		private T a4;
		private T a6;
		private T b2;
		private T b4;
		private T b6;
		private T b8;
		private T c4;
		private T c6;
		private T discriminant;
		private T j;
	}

	private static <T extends Element<T>> CurveParameters<T> computeCurveParameters(Field<T> field, T a, T b) {
		return computeCurveParameters(field, field.zero(), field.zero(), field.zero(), a, b);
	}

	// y^2 + a1 xy + a3 y = x^3 + a2 x^2 + a4 x + a6
	private static <T extends Element<T>> CurveParameters<T> computeCurveParameters(Field<T> field, T a1, T a2, T a3,
			T a4, T a6) {
		CurveParameters<T> result = new CurveParameters<>();
		result.a1 = a1;
		result.a2 = a2;
		result.a3 = a3;
		result.a4 = a4;
		result.a6 = a6;
		result.b2 = field.add(field.multiply(a1, a1), field.multiply(4, a2));
		result.b4 = field.add(field.multiply(a1, a3), field.multiply(2, a4));
		result.b6 = field.add(field.multiply(a3, a3), field.multiply(4, a6));
		result.b8 = field.add(field.add(field.multiply(a1, a1, a6), field.multiply(-1, a1, a3, a4)),
				field.multiply(4, a2, a6), field.multiply(a2, a3, a3), field.multiply(-1, a4, a4));
		result.c4 = field.add(field.multiply(result.b2, result.b2), field.multiply(-24, result.b4));
		result.c6 = field.add(field.multiply(-1, result.b2, result.b2, result.b2),
				field.multiply(36, result.b2, result.b4), field.multiply(-216, result.b6));
		result.discriminant = field.add(field.multiply(-1, result.b2, result.b2, result.b8),
				field.multiply(-8, result.b4, result.b4, result.b4), field.multiply(-27, result.b6, result.b6),
				field.multiply(9, result.b2, result.b4, result.b6));
		result.j = field.divide(field.power(result.c4, 3), result.discriminant);
		return result;
	}

	private static <T extends Element<T>> Polynomial<T> polynomialFromParameters(Field<T> field,
			CurveParameters<T> parameters) {
		PolynomialRing<T> r = AbstractPolynomialRing.getPolynomialRing(field, 2, Monomial.REVLEX);
		Polynomial<T> result = r.add(r.getVarPower(1, 3), r.getEmbedding(parameters.a2, new int[] { 2, 0 }),
				r.getEmbedding(parameters.a4, new int[] { 1, 0 }), r.getEmbedding(parameters.a6));
		result = r.subtract(result, r.add(r.getVarPower(2, 2), r.getEmbedding(parameters.a1, new int[] { 1, 1 }),
				r.getEmbedding(parameters.a3, new int[] { 0, 1 })));
		return result;
	}

	private static <T extends Element<T>> CurveParameters<T> parametersFromPolynomial(Field<T> field,
			Polynomial<T> polynomial) {
		PolynomialRing<T> r = AbstractPolynomialRing.getPolynomialRing(field, 2, Monomial.REVLEX);
		T a1 = field.negative(polynomial.coefficient(r.getMonomial(new int[] { 1, 1 })));
		T a2 = polynomial.coefficient(r.getMonomial(new int[] { 2, 0 }));
		T a3 = field.negative(polynomial.coefficient(r.getMonomial(new int[] { 1, 1 })));
		;
		T a4 = polynomial.coefficient(r.getMonomial(new int[] { 1, 0 }));
		;
		T a6 = polynomial.coefficient(r.getMonomial(new int[] { 0, 0 }));
		;
		return computeCurveParameters(field, a1, a2, a3, a4, a6);
	}

	// x' = x - xOffset, y' = y - x * yxOffset - yOffset
	private static <T extends Element<T>> CurveParameters<T> substituteVariables(Field<T> field,
			CurveParameters<T> original, T xOffset, T yOffset, T yxOffset) {
		PolynomialRing<T> r = AbstractPolynomialRing.getPolynomialRing(field, 2, Monomial.REVLEX);
		Polynomial<T> asPolynomial = polynomialFromParameters(field, original);
		List<Polynomial<T>> substitute = new ArrayList<>();
		substitute.add(r.add(r.getVar(1), r.getEmbedding(xOffset)));
		substitute.add(r.add(r.getVar(2), r.getEmbedding(yxOffset, new int[] { 1, 0 }), r.getEmbedding(yOffset)));
		Polynomial<T> substituted = r.substitute(asPolynomial, substitute);
		return parametersFromPolynomial(field, substituted);
	}

	public static class ReductionResult<T extends Element<T>, S extends Element<S>> {
		private EllipticCurve<T> curve;
		private KodairaSymbol kodairaSymbol;
		private int conductorExponent;
		private int localIndex;
		private EllipticCurve<S> reducedCurve;
		private MathMap<ProjectivePoint<T>, ProjectivePoint<S>> reductionMap;
		private MathMap<ProjectivePoint<S>, ProjectivePoint<T>> liftAsTorsionPointMap;

		private ReductionResult(EllipticCurve<T> curve, EllipticCurve<S> reducedCurve,
				MathMap<ProjectivePoint<T>, ProjectivePoint<S>> reductionMap,
				MathMap<ProjectivePoint<S>, ProjectivePoint<T>> liftAsTorsionPointMap) {
			this.curve = curve;
			this.reducedCurve = reducedCurve;
			this.reductionMap = reductionMap;
			this.liftAsTorsionPointMap = liftAsTorsionPointMap;
		}

		public EllipticCurve<T> getCurve() {
			return curve;
		}

		public EllipticCurve<S> getReducedCurve() {
			return reducedCurve;
		}

		public MathMap<ProjectivePoint<T>, ProjectivePoint<S>> getReductionMap() {
			return reductionMap;
		}

		public MathMap<ProjectivePoint<S>, ProjectivePoint<T>> getLiftAsTorsionPointMap() {
			return liftAsTorsionPointMap;
		}
	}

	public static <T extends Element<T>, S extends Element<S>> ReductionResult<T, S> reduce(
			DiscreteValuationField<T, S> field, EllipticCurve<T> curve) {
		Field<S> reduction = field.residueField();
		UnivariatePolynomialRing<S> reductionPolynomials = reduction.getUnivariatePolynomialRing();
		CurveParameters<T> parameters = computeCurveParameters(field, curve.getA(), curve.getB());
		Value discriminantValue = field.valuation(parameters.discriminant);
		if (discriminantValue.equals(Value.ZERO)) {
			EllipticCurve<S> reducedCurve = new EllipticCurve<S>(reduction, field.reduceInteger(curve.getA()),
					field.reduceInteger(curve.getB()));
			KodairaSymbol symbol = KodairaSymbol.I0;
			int conductorExponent = 0;
			int localIndex = 1;
		}
		S reducedA = field.reduceInteger(curve.getA());
		S reducedB = field.reduceInteger(curve.getB());
		List<S> reductionCoefficients = new ArrayList<>();
		reductionCoefficients.add(reducedB);
		reductionCoefficients.add(reducedA);
		reductionCoefficients.add(reduction.zero());
		reductionCoefficients.add(reduction.one());
		UnivariatePolynomial<S> reducedRhs = reductionPolynomials.getPolynomial(reductionCoefficients);
		FactorizationResult<Polynomial<S>, S> squareFreeFactors = reductionPolynomials
				.squareFreeFactorization(reducedRhs);
		for (Polynomial<S> squareFreeFactor : squareFreeFactors.primeFactors()) {
			if (squareFreeFactors.multiplicity(squareFreeFactor) > 1) {
				S xOffsetReduced = reduction
						.negative(reductionPolynomials.toUnivariate(squareFreeFactor).univariateCoefficient(0));
				parameters = substituteVariables(field, parameters, field.liftToInteger(xOffsetReduced), field.zero(),
						field.zero());
			}
		}
		if (field.valuation(parameters.b2).equals(Value.ZERO)) {
			KodairaSymbol symbol = KodairaSymbol.Iv;
			int conductorExponent = 1;
			int localIndex = 1;
			int value = field.valuation(parameters.discriminant).value();
		}
		if (field.valuation(parameters.a6).compareTo(new Value(2)) < 0) {
			KodairaSymbol symbol = KodairaSymbol.II;
			int conductorExponent = 1;
			int localIndex = field.valuation(parameters.discriminant).value();
		}
		if (field.valuation(parameters.b8).compareTo(new Value(3)) < 3) {
			KodairaSymbol symbol = KodairaSymbol.III;
			int conductorExponent = 2;
			int localIndex = field.valuation(parameters.discriminant).value() - 1;

		}
		if (field.valuation(parameters.b6).compareTo(new Value(3)) < 3) {
			KodairaSymbol symbol = KodairaSymbol.IV;
			int conductorExponent = 3;
			int localIndex = field.valuation(parameters.discriminant).value() - 2;

		}
		EllipticCurve<S> reducedCurve = new EllipticCurve<S>(reduction, field.reduceInteger(curve.getA()),
				field.reduceInteger(curve.getB()));
		MathMap<ProjectivePoint<T>, ProjectivePoint<S>> reductionMap = new MathMap<>() {

			@Override
			public ProjectivePoint<S> evaluate(ProjectivePoint<T> t) {
				List<S> result = new ArrayList<>();
				for (T coordinate : t.getCoords()) {
					result.add(field.reduceInteger(coordinate));
				}
				return new ProjectivePoint<>(reduction, result);
			}
		};
		MathMap<ProjectivePoint<S>, ProjectivePoint<T>> liftAsTorsionPointMap = new MathMap<>() {

			@Override
			public ProjectivePoint<T> evaluate(ProjectivePoint<S> t) {
				int order = reducedCurve.getOrder(t).intValueExact();
				if (order == 1) {
					return curve.neutral();
				}
				if (order == 2) {
					T x = field.ringOfIntegers().henselLift(field.getUnivariatePolynomialRing().toUnivariate(curve.rhs),
							t.getDehomogenisedCoord(1, 3));
					return new ProjectivePoint<>(field, x, field.zero(), field.one());
				}
				Polynomial<T> divisionPolynomial = curve.getDivisionPolynomial(order);
				if (order % 2 == 0) {
					divisionPolynomial = curve.affineRing.divide(divisionPolynomial, curve.affineRing.getVar(2));
				}
				UnivariatePolynomial<T> xOnlyDivisionPolynomial = field.getUnivariatePolynomialRing()
						.getEmbedding(divisionPolynomial, new int[] { 0 });
				T x = field.ringOfIntegers().henselLift(xOnlyDivisionPolynomial, t.getDehomogenisedCoord(1, 3));
				List<T> eval = new ArrayList<>();
				eval.add(x);
				eval.add(null);
				UnivariatePolynomial<T> yPolynomial = field.getUnivariatePolynomialRing().getEmbedding(
						curve.affineRing.partiallyEvaluate(curve.definingPolynomial, eval), new int[] { -1, 0 });
				T y = field.ringOfIntegers().henselLift(yPolynomial, t.getDehomogenisedCoord(2, 3));
				return new ProjectivePoint<>(field, x, y, field.one());
			}
		};
		return new ReductionResult<>(curve, reducedCurve, reductionMap, liftAsTorsionPointMap);
	}

	public static IntE numberOfSupersingularCurves(IntE prime) {
		Integers z = Integers.z();
		QuotientAndRemainderResult<IntE> qr12 = z.quotientAndRemainder(prime, z.getInteger(12));
		IntE add = null;
		if (qr12.getRemainder().equals(z.one())) {
			add = z.zero();
		} else if (qr12.getRemainder().equals(z.getInteger(5)) || qr12.getRemainder().equals(z.getInteger(7))) {
			add = z.one();
		} else if (qr12.getRemainder().equals(z.getInteger(11))) {
			add = z.getInteger(2);
		}
		return z.add(qr12.getQuotient(), add);
	}

	public EllipticCurve(Field<T> field, T a, T b) {
		this(field, field.zero(), field.zero(), field.zero(), a, b);
	}

	public EllipticCurve(Field<T> field, T a1, T a2, T a3, T a4, T a6) {
		super(asGenericProjectiveScheme(field, a1, a2, a3, a4, a6));
		this.field = field;
		this.a1 = a1;
		this.a2 = a2;
		this.a3 = a3;
		this.a4 = a4;
		this.a6 = a6;
		this.b2 = field.add(field.multiply(a1, a1), field.multiply(4, a2));
		this.b4 = field.add(field.multiply(a1, a3), field.multiply(2, a4));
		this.b6 = field.add(field.multiply(a3, a3), field.multiply(4, a6));
		this.b8 = field.add(field.add(field.multiply(a1, a1, a6), field.multiply(-1, a1, a3, a4)),
				field.multiply(4, a2, a6), field.multiply(a2, a3, a3), field.multiply(-1, a4, a4));
		this.c4 = field.add(field.multiply(b2, b2), field.multiply(-24, b4));
		this.c6 = field.add(field.multiply(-1, b2, b2, b2), field.multiply(36, b2, b4), field.multiply(-216, b6));
		this.discriminant = field.add(field.multiply(-1, b2, b2, b8), field.multiply(-8, b4, b4, b4),
				field.multiply(-27, b6, b6), field.multiply(9, b2, b4, b6));
		if (discriminant.equals(this.field.zero())) {
			throw new ArithmeticException("Singular curve");
		}
		this.j = field.divide(field.power(c4, 3), discriminant);
		this.pointAtInfinity = new ProjectivePoint<T>(this.field, this.field.zero(), this.field.one(),
				this.field.zero());
		this.projectiveRing = asGenericProjectiveScheme().homogenousPolynomialRing();
		this.affineRing = AbstractPolynomialRing.getPolynomialRing(this.field, 2, Monomial.REVLEX);
		PolynomialRing<T> r = this.affineRing;
		Polynomial<T> a6p = r.getEmbedding(a6);
		Polynomial<T> a4x = r.getEmbedding(a4, new int[] { 1, 0 });
		Polynomial<T> a2x2 = r.getEmbedding(a2, new int[] { 2, 0 });
		Polynomial<T> x3 = r.getVarPower(1, 3);
		Polynomial<T> y2 = r.getVarPower(2, 2);
		Polynomial<T> a1xy = r.getEmbedding(a1, new int[] { 1, 1 });
		Polynomial<T> a3y = r.getEmbedding(a3, new int[] { 0, 1 });
		this.rhs = r.add(x3, a2x2, a4x, a6p);
		this.lhs = r.add(y2, a1xy, a3y);
		this.definingPolynomial = r.subtract(rhs, lhs);
		this.univariateRing = this.field.getUnivariatePolynomialRing();
		this.torsionPointBasis = new TreeMap<>();
	}

	public T jInvariant() {
		return j;
	}

	public IdentityIsomorphism<T> identity() {
		return new IdentityIsomorphism<>(this);
	}

	private List<Polynomial<T>> linearTransforms(PolynomialRing<T> r) {
		List<Polynomial<T>> equations = new ArrayList<>();
		equations.add(r.add(r.multiply(2, r.getVarPower(1, 3), r.getVar(3)), r.multiply(a1, r.getVar(1))));
		equations.add(r.add(r.multiply(-1, r.getVarPower(1, 6), r.getVarPower(3, 2)),
				r.multiply(-1, r.getEmbedding(a1), r.getVarPower(1, 4), r.getVar(3)),
				r.multiply(a2, r.getVarPower(1, 2)), r.multiply(3, r.getVarPower(1, 2), r.getVar(2))));
		equations.add(r.add(r.multiply(r.getEmbedding(a1), r.getVarPower(1, 3), r.getVar(2)),
				r.multiply(a3, r.getVarPower(1, 3)), r.multiply(2, r.getVarPower(1, 3), r.getVar(4))));
		equations.add(r.add(
				r.add(r.multiply(r.multiply(-1, r.getEmbedding(a1), r.getVarPower(1, 6)), r.getVar(2), r.getVar(3)),
						r.multiply(-1, r.getEmbedding(a3), r.multiply(r.getVarPower(1, 6), r.getVar(3)))),
				r.add(r.multiply(-2, r.getVarPower(1, 6), r.getVar(3), r.getVar(4)),
						r.multiply(2, r.getEmbedding(a2), r.getVarPower(1, 4), r.getVar(2))),
				r.add(r.multiply(3, r.getVarPower(1, 4), r.getVarPower(2, 2)),
						r.multiply(-1, r.getEmbedding(a1), r.getVarPower(1, 4), r.getVar(4)),
						r.multiply(a4, r.getVarPower(1, 4)))));
		equations.add(r.add(
				r.add(r.multiply(r.getEmbedding(a2), r.getVarPower(1, 6), r.getVarPower(2, 2)),
						r.multiply(r.getVarPower(1, 6), r.getVarPower(2, 3)),
						r.multiply(r.multiply(-1, r.getEmbedding(a1), r.getVarPower(1, 6)), r.getVar(2), r.getVar(4))),
				r.add(r.multiply(r.getEmbedding(a4), r.getVarPower(1, 6), r.getVar(2)),
						r.multiply(-1, r.getEmbedding(a3), r.getVarPower(1, 6), r.getVar(4))),
				r.add(r.multiply(-1, r.getVarPower(1, 6), r.getVarPower(4, 2)), r.multiply(a6, r.getVarPower(1, 6)))));
		return equations;
	}

	public List<Isogeny<T>> getIsomorphisms(EllipticCurve<T> range) {
		if (!this.jInvariant().equals(range.jInvariant())) {
			throw new ArithmeticException("Not isomorphic");
		}
//		a1' = 2*s^3*xyO + a1*s
//		a2' = -1*s^6*xyO^2 + -1*a1*s^4*xyO + a2*s^2 + 3*s^2*xO
//		a3' = a1*s^3*xO + a3*s^3 + 2*s^3yO
//		a4' = -1*a1s^6*xO*xyO + -1*a3*s^6*xyO + -2*s^6*xyO*yO + 2*a2*s^4*xO + 3*s^4*xO^2 + -1*a1*s^4*yO + a4*s^4
//		a6' = a2*s^6*xO^2 + s^6*xO^3 + -1*a1*s^6*xO*yO + a4*s^6*xO + -1*a3*s^6*yO + -1*s^6*yO^2 + a6*s^6
		PolynomialRing<T> r = AbstractPolynomialRing.getPolynomialRing(field, 4, Monomial.REVLEX);
		List<Polynomial<T>> equations = linearTransforms(r);
		equations.set(0, r.subtract(equations.get(0), r.getEmbedding(a1)));
		equations.set(1, r.subtract(equations.get(1), r.getEmbedding(a2)));
		equations.set(2, r.subtract(equations.get(2), r.getEmbedding(a3)));
		equations.set(3, r.subtract(equations.get(3), r.getEmbedding(a4)));
		equations.set(4, r.subtract(equations.get(4), r.getEmbedding(a6)));
		List<AffinePoint<T>> solved = r.solve(equations);
		List<Isogeny<T>> isomorphisms = new ArrayList<>();
		for (AffinePoint<T> root : solved) {
			isomorphisms.add(new Isomorphism<>(this, field.inverse(root.getCoord(1)), root.getCoord(2),
					root.getCoord(3), root.getCoord(4)));
		}
		return isomorphisms;
	}

	public List<Isogeny<T>> getAutomorphisms() {
		return getIsomorphisms(this);
	}

	public Isogeny<T> frobeniousEndomorphism() {
		T j = jInvariant();
		BigInteger p = field.characteristic();
		BigInteger q = field.getNumberOfElements();
		int n = 0;
		T jp = j;
		while (!q.equals(BigInteger.ONE)) {
			n++;
			q = q.divide(p);
			jp = field.power(jp, p);
			if (j.equals(jp)) {
				break;
			}
		}
		Isogeny<T> frobenius = new Frobenius<>(this, n);
		if (frobenius.getRange().equals(this)) {
			return frobenius;
		}
		List<Isogeny<T>> isomorphisms = frobenius.getRange().getIsomorphisms(this);
		return new CompositionIsogeny<>(frobenius, isomorphisms.get(0));
	}

	public Graph<EllipticCurve<T>, Isogeny<T>> isogenyGraph(int prime) {
		if (!getField().isFinite()) {
			throw new UnsupportedOperationException("Not implemented");
		}
		GraphBuilder<EllipticCurve<T>, Isogeny<T>> builder = new GraphBuilder<>(true);
		Set<T> visited = new TreeSet<>();
		Map<T, EllipticCurve<T>> canonical = new TreeMap<>();
		Deque<EllipticCurve<T>> queue = new LinkedList<>();
		queue.add(this);
		while (!queue.isEmpty()) {
			EllipticCurve<T> curve = queue.poll();
			if (visited.contains(curve.jInvariant())) {
				continue;
			}
			builder.addVertex(curve);
			visited.add(curve.jInvariant());
			canonical.put(curve.jInvariant(), curve);
			List<ProjectivePoint<T>> torsionPoints = curve.getTorsionPointBasis(prime);
			List<ProjectivePoint<T>> kernelPoints = new ArrayList<>();
			if (torsionPoints.isEmpty()) {
				continue;
			}
			kernelPoints.add(torsionPoints.get(0));
			if (torsionPoints.size() == 2) {
				for (int i = 0; i < prime; i++) {
					kernelPoints.add(curve.add(curve.multiply(i, torsionPoints.get(0)), torsionPoints.get(1)));
				}
			}
			for (ProjectivePoint<T> torsionPoint : kernelPoints) {
				Isogeny<T> isogeny = new KernelPointIsogeny<>(curve, torsionPoint, prime);
				EllipticCurve<T> range = isogeny.getRange();
				if (canonical.containsKey(range.jInvariant())) {
					Isogeny<T> isomorphism = range.getIsomorphisms(canonical.get(range.jInvariant())).get(0);
					isogeny = new CompositionIsogeny<>(isogeny, isomorphism);
					range = canonical.get(range.jInvariant());
				} else {
					canonical.put(range.jInvariant(), range);
				}
				builder.addEdge(curve, range, isogeny);
				if (!visited.contains(range.jInvariant())) {
					queue.add(range);
				}
			}
		}
		return builder.build();

	}

	public BigInteger trace() {
		return field.getNumberOfElements().add(BigInteger.ONE).subtract(getNumberOfElements());
	}

	public BigInteger asMultiplicationIsogenyModOrder(Isogeny<T> isogeny) {
		Integers z = Integers.z();
		BigInteger numberOfPoints = getNumberOfElements();
		FactorizationResult<IntE, IntE> factorization = z.uniqueFactorization(z.getInteger(numberOfPoints));
		List<IntE> primePowers = new ArrayList<>();
		List<IntE> residues = new ArrayList<>();
		for (IntE prime : factorization.primeFactors()) {
			int exp = MiscAlgorithms.DivRoundUp(factorization.multiplicity(prime), 2);
			IntE power = z.power(prime, exp);
			IntE subPower = z.divideChecked(power, prime);
			IntE cofactor = z.divideChecked(z.getInteger(numberOfPoints),
					z.power(prime, factorization.multiplicity(prime)));
			ProjectivePoint<T> torsionPoint;
			do {
				torsionPoint = multiply(cofactor, getRandomElement());
			} while (multiply(subPower, torsionPoint).equals(neutral()));
			while (!multiply(power, torsionPoint).equals(neutral())) {
				power = z.multiply(prime, power);
				exp++;
			}
			primePowers.add(power);
			IntE currentPower = power;
			IntE currentReconstructed = z.zero();
			IntE currentDigit = z.one();
			do {
				currentPower = z.divideChecked(power, prime);
				currentDigit = z.multiply(currentDigit, prime);
				ProjectivePoint<T> point = multiply(currentPower, torsionPoint);
				ProjectivePoint<T> image = isogeny.evaluate(point);
				for (int i = 0; i < prime.intValueExact(); i++) {
					ProjectivePoint<T> multiplied = multiply(i, point);
					if (multiplied.equals(image)) {
						currentReconstructed = z.add(z.multiply(z.getInteger(i), currentDigit), currentReconstructed);
						break;
					}
				}
			} while (!currentPower.equals(z.one()));
			residues.add(currentReconstructed);
		}
		return z.chineseRemainderTheoremModuli(residues, primePowers).getValue();
	}

	public boolean isSupersingular() {
		if (!superSingularDeterminied) {
			this.superSingular = attemptSupersingularCount();
			this.superSingularDeterminied = true;
		}
		return superSingular;
	}

	public T getA1() {
		return a1;
	}

	public T getA2() {
		return a2;
	}

	public T getA3() {
		return a3;
	}

	public T getA4() {
		return a4;
	}

	public T getA6() {
		return a6;
	}

	public boolean isWeierstrass() {
		if (field.characteristic().equals(BigInteger.TWO)) {
			if (a1.equals(field.zero())) {
				return a2.equals(field.zero());
			} else {
				return a1.equals(field.one()) && a3.equals(field.zero()) && a4.equals(field.zero());
			}
		}
		if (!a1.equals(field.zero()) || !a3.equals(field.zero()))
			return false;
		if (field.characteristic().equals(BigInteger.valueOf(3))) {
			return true;
		}
		return a2.equals(field.zero());
	}

	public Isogeny<T> getWeierstrassForm() {
		if (weierstrassForm != null) {
			return weierstrassForm;
		}
		T zero = field.zero();
		T one = field.one();
		if (field.characteristic().equals(BigInteger.TWO)) {
			if (!a1.equals(zero)) {
				if (a1.equals(one) && a3.equals(zero) && a4.equals(zero)) {
					weierstrassForm = new IdentityIsomorphism<>(this);
					return weierstrassForm;
				}
				PolynomialRing<T> r4 = AbstractPolynomialRing.getPolynomialRing(field, 4, Monomial.REVLEX);
				PolynomialRing<T> r3 = AbstractPolynomialRing.getPolynomialRing(field, 3, Monomial.REVLEX);
				List<T> eval = new ArrayList<>();
				eval.add(null);
				eval.add(null);
				eval.add(zero);
				eval.add(null);
				int[] map = new int[] { 0, 1, -1, 2 };
				List<Polynomial<T>> equations = linearTransforms(r4);
				List<Polynomial<T>> requiredEquations = new ArrayList<>();
				requiredEquations
						.add(r3.subtract(r3.getEmbedding(r4.partiallyEvaluate(equations.get(0), eval), map), r3.one()));
				requiredEquations.add(r3.getEmbedding(r4.partiallyEvaluate(equations.get(3), eval), map));
				requiredEquations.add(r3.getEmbedding(r4.partiallyEvaluate(equations.get(2), eval), map));
				List<AffinePoint<T>> transform = r3.solve(requiredEquations);
				weierstrassForm = new Isomorphism<>(this, field.inverse(transform.get(0).getCoord(1)),
						transform.get(0).getCoord(2), zero, transform.get(0).getCoord(3));
				return weierstrassForm;
			} else {
				PolynomialRing<T> r4 = AbstractPolynomialRing.getPolynomialRing(field, 4, Monomial.REVLEX);
				PolynomialRing<T> r1 = AbstractPolynomialRing.getPolynomialRing(field, 1, Monomial.REVLEX);
				List<T> eval = new ArrayList<>();
				eval.add(one);
				eval.add(null);
				eval.add(zero);
				eval.add(zero);
				int[] map = new int[] { -1, 0, -1, -1 };
				List<Polynomial<T>> equations = linearTransforms(r4);
				List<Polynomial<T>> requiredEquations = new ArrayList<>();
				requiredEquations.add(r1.getEmbedding(r4.partiallyEvaluate(equations.get(1), eval), map));
				List<AffinePoint<T>> transform = r1.solve(requiredEquations);
				weierstrassForm = new Isomorphism<>(this, one, transform.get(0).getCoord(1), zero, zero);
				return weierstrassForm;
			}
		} else if (field.characteristic().equals(BigInteger.valueOf(3))) {
			if (a1.equals(zero) && a3.equals(zero)) {
				weierstrassForm = new IdentityIsomorphism<>(this);
				return weierstrassForm;
			}
			PolynomialRing<T> r4 = AbstractPolynomialRing.getPolynomialRing(field, 4, Monomial.REVLEX);
			PolynomialRing<T> r2 = AbstractPolynomialRing.getPolynomialRing(field, 2, Monomial.REVLEX);
			List<T> eval = new ArrayList<>();
			eval.add(one);
			eval.add(zero);
			eval.add(null);
			eval.add(null);
			int[] map = new int[] { -1, -1, 0, 1 };
			List<Polynomial<T>> equations = linearTransforms(r4);
			List<Polynomial<T>> requiredEquations = new ArrayList<>();
			requiredEquations.add(r2.getEmbedding(r4.partiallyEvaluate(equations.get(0), eval), map));
			requiredEquations.add(r2.getEmbedding(r4.partiallyEvaluate(equations.get(2), eval), map));
			List<AffinePoint<T>> transform = r2.solve(requiredEquations);
			weierstrassForm = new Isomorphism<>(this, one, zero, transform.get(0).getCoord(1),
					transform.get(0).getCoord(2));
			return weierstrassForm;
		}
		if (a1.equals(zero) && a3.equals(zero) && a2.equals(zero)) {
			weierstrassForm = new IdentityIsomorphism<>(this);
			return weierstrassForm;
		}
		PolynomialRing<T> r4 = AbstractPolynomialRing.getPolynomialRing(field, 4, Monomial.REVLEX);
		PolynomialRing<T> r3 = AbstractPolynomialRing.getPolynomialRing(field, 3, Monomial.REVLEX);
		List<T> eval = new ArrayList<>();
		eval.add(one);
		eval.add(null);
		eval.add(null);
		eval.add(null);
		int[] map = new int[] { -1, 0, 1, 2 };
		List<Polynomial<T>> equations = linearTransforms(r4);
		List<Polynomial<T>> requiredEquations = new ArrayList<>();
		requiredEquations.add(r3.getEmbedding(r4.partiallyEvaluate(equations.get(0), eval), map));
		requiredEquations.add(r3.getEmbedding(r4.partiallyEvaluate(equations.get(1), eval), map));
		requiredEquations.add(r3.getEmbedding(r4.partiallyEvaluate(equations.get(2), eval), map));
		List<AffinePoint<T>> transform = r3.solve(requiredEquations);
		weierstrassForm = new Isomorphism<>(this, one, transform.get(0).getCoord(1), transform.get(0).getCoord(2),
				transform.get(0).getCoord(3));
		return weierstrassForm;
	}

	public EllipticCurve<T> getQuadraticTwist() {
		if (field.characteristic().equals(BigInteger.TWO)) {
			T twist;
			UnivariatePolynomial<T> twistPolynomial;
			do {
				twist = field.getRandomElement();
				List<T> coeffs = new ArrayList<>();
				coeffs.add(twist);
				coeffs.add(field.one());
				coeffs.add(field.one());
				twistPolynomial = univariateRing.getPolynomial(coeffs);
			} while (field.hasRoots(twistPolynomial));
			return new EllipticCurve<>(field, a1, field.add(a2, field.multiply(twist, a1, a1)), a3, a4,
					field.add(a6, field.multiply(twist, a3, a3)));
		}
		if (!a1.equals(field.zero()) || !a3.equals(field.zero())) {
			PolynomialRing<T> r4 = AbstractPolynomialRing.getPolynomialRing(field, 4, Monomial.REVLEX);
			PolynomialRing<T> r2 = AbstractPolynomialRing.getPolynomialRing(field, 2, Monomial.REVLEX);
			List<T> eval = new ArrayList<>();
			eval.add(field.one());
			eval.add(field.zero());
			eval.add(null);
			eval.add(null);
			int[] map = new int[] { -1, -1, 0, 1 };
			List<Polynomial<T>> equations = linearTransforms(r4);
			List<Polynomial<T>> requiredEquations = new ArrayList<>();
			requiredEquations.add(r2.getEmbedding(r4.partiallyEvaluate(equations.get(0), eval), map));
			requiredEquations.add(r2.getEmbedding(r4.partiallyEvaluate(equations.get(2), eval), map));
			List<AffinePoint<T>> transform = r2.solve(requiredEquations);
			return new Isomorphism<>(this, field.one(), field.zero(), transform.get(0).getCoord(1),
					transform.get(0).getCoord(2)).getRange().getQuadraticTwist();

		}
		T twist;
		do {
			twist = field.getRandomElement();
		} while (field.hasSqrt(twist));
		return new EllipticCurve<>(field, field.zero(), field.multiply(twist, a2), field.zero(),
				field.multiply(twist, twist, a4), field.multiply(field.power(twist, 3), a6));
	}

	private T finiteFieldTrace(T t) {
		return finiteFieldTrace(t, 1);
	}

	private T finiteFieldTrace(T t, int degree) {
		BigInteger q = field.getNumberOfElements();
		BigInteger p = field.characteristic().pow(degree);
		T result = field.zero();
		T power = t;
		while (!q.equals(BigInteger.ONE)) {
			result = field.add(power, result);
			power = field.power(power, p);
			q = q.divide(p);
		}
		return result;
	}

	private T finiteFieldNorm(T t) {
		return finiteFieldNorm(t, 1);
	}

	private T finiteFieldNorm(T t, int degree) {
		BigInteger q = field.getNumberOfElements();
		BigInteger p = field.characteristic().pow(degree);
		T result = field.one();
		T power = t;
		while (!q.equals(BigInteger.ONE)) {
			result = field.multiply(power, result);
			power = field.power(power, p);
			q = q.divide(p);
		}
		return result;
	}

	public T getA() {
		if (field.characteristic().equals(BigInteger.valueOf(2))
				|| field.characteristic().equals(BigInteger.valueOf(3))) {
			throw new ArithmeticException("Not defined!");
		}
		return getWeierstrassForm().getRange().getA4();
	}

	public T getB() {
		if (field.characteristic().equals(BigInteger.valueOf(2))
				|| field.characteristic().equals(BigInteger.valueOf(3))) {
			throw new ArithmeticException("Not defined!");
		}
		return getWeierstrassForm().getRange().getA6();
	}

	public CoordinateRing<T> getCoordinateRing() {
		if (this.coordinateRing == null) {
			this.coordinateRing = affineRing.getIdeal(Collections.singletonList(definingPolynomial)).divideOut();
		}
		return this.coordinateRing;
	}

	@Override
	public boolean hasRationalPoint(ProjectivePoint<T> p) {
		if (p.getDim() != 2)
			return false;
		if (p.equals(this.pointAtInfinity))
			return true;
		if (p.getCoord(3).equals(this.field.zero()))
			return false;
		return this.affineRing.evaluate(this.definingPolynomial, p.getDehomogenous(3).getCoords())
				.equals(this.field.zero());
	}

	public Polynomial<T> getDivisionPolynomial(int m) {
		if (this.divisionPolynomials == null) {
			this.divisionPolynomials = new TreeMap<>();
		}
		if (!this.divisionPolynomials.containsKey(m)) {
			CoordinateRing<T> r = getCoordinateRing();
			if (m == 0) {
				divisionPolynomials.put(m, r.zero().getElement());
			} else if (m == 1) {
				divisionPolynomials.put(m, r.one().getElement());
			} else if (m == 2) {
				CoordinateRingElement<T> twoY = r.multiply(2, r.getVar(2));
				CoordinateRingElement<T> a1x = r.multiply(a1, r.getVar(1));
				CoordinateRingElement<T> a3cr = r.getEmbedding(affineRing.getEmbedding(a3));
				divisionPolynomials.put(m, r.add(twoY, a1x, a3cr).getElement());
			} else if (m == 3) {
				CoordinateRingElement<T> threeX4 = r.multiply(3, r.power(r.getVar(1), 4));
				CoordinateRingElement<T> b2x3 = r.multiply(b2, r.power(r.getVar(1), 3));
				CoordinateRingElement<T> threeB4X2 = r.multiply(3, r.getEmbedding(affineRing.getEmbedding(b4)),
						r.power(r.getVar(1), 2));
				CoordinateRingElement<T> threeB6X = r.multiply(3, r.getEmbedding(affineRing.getEmbedding(b6)),
						r.getVar(1));
				CoordinateRingElement<T> b8cr = r.getEmbedding(affineRing.getEmbedding(b8));
				divisionPolynomials.put(m, r.add(r.add(threeX4, b2x3, threeB4X2), r.add(threeB6X, b8cr)).getElement());
			} else if (m == 4) {
				CoordinateRingElement<T> psi2 = r.getEmbedding(getDivisionPolynomial(2));
				CoordinateRingElement<T> twoX6 = r.multiply(2, r.power(r.getVar(1), 6));
				CoordinateRingElement<T> b2x5 = r.multiply(b2, r.power(r.getVar(1), 5));
				CoordinateRingElement<T> fiveB4X4 = r.multiply(5, r.getEmbedding(affineRing.getEmbedding(b4)),
						r.power(r.getVar(1), 4));
				CoordinateRingElement<T> tenB6X3 = r.multiply(10, r.getEmbedding(affineRing.getEmbedding(b6)),
						r.power(r.getVar(1), 3));
				CoordinateRingElement<T> tenB8X2 = r.multiply(10, r.getEmbedding(affineRing.getEmbedding(b8)),
						r.power(r.getVar(1), 2));
				CoordinateRingElement<T> b2b8mb4b6x = r
						.multiply(field.subtract(field.multiply(b2, b8), field.multiply(b4, b6)), r.getVar(1));
				CoordinateRingElement<T> b4b8mb6sqr = r.getEmbedding(
						affineRing.getEmbedding(field.subtract(field.multiply(b4, b8), field.multiply(b6, b6))));
				divisionPolynomials.put(m,
						r.multiply(psi2,
								r.add(r.add(twoX6, b2x5, fiveB4X4, tenB6X3), r.add(tenB8X2, b2b8mb4b6x, b4b8mb6sqr)))
								.getElement());
			} else if (m % 2 == 0) {
				int n = m / 2;
				CoordinateRingElement<T> psiNm2 = r.getEmbedding(getDivisionPolynomial(n - 2));
				CoordinateRingElement<T> psiNm1 = r.getEmbedding(getDivisionPolynomial(n - 1));
				CoordinateRingElement<T> psiN = r.getEmbedding(getDivisionPolynomial(n));
				CoordinateRingElement<T> psiNp1 = r.getEmbedding(getDivisionPolynomial(n + 1));
				CoordinateRingElement<T> psiNp2 = r.getEmbedding(getDivisionPolynomial(n + 2));
				CoordinateRingElement<T> psiMpsi2 = r.subtract(r.multiply(r.power(psiNm1, 2), psiN, psiNp2),
						r.multiply(psiNm2, psiN, r.power(psiNp1, 2)));
				Polynomial<T> psiM = r.divideChecked(psiMpsi2, r.getEmbedding(getDivisionPolynomial(2))).getElement();
				divisionPolynomials.put(m, psiM);
			} else {
				int n = (m - 1) / 2;
				CoordinateRingElement<T> psiNm1 = r.getEmbedding(getDivisionPolynomial(n - 1));
				CoordinateRingElement<T> psiN = r.getEmbedding(getDivisionPolynomial(n));
				CoordinateRingElement<T> psiNp1 = r.getEmbedding(getDivisionPolynomial(n + 1));
				CoordinateRingElement<T> psiNp2 = r.getEmbedding(getDivisionPolynomial(n + 2));
				divisionPolynomials.put(m,
						r.subtract(r.multiply(psiNp2, r.power(psiN, 3)), r.multiply(psiNm1, r.power(psiNp1, 3)))
								.getElement());
			}
		}
		return this.divisionPolynomials.get(m);
//			return p;
//			if (m % 2 == 0) {
//				return this.affineRing.multiply(2, this.affineRing.getVar(2), p);
//			} else {
//				return p;
//			}
//		}
//		PolynomialRing<T> r = this.affineRing;
//		Polynomial<T> o = r.one();
//		Polynomial<T> x = r.getVar(1);
//		Polynomial<T> x2 = r.getVarPower(1, 2);
//		Polynomial<T> x3 = r.getVarPower(1, 3);
//		Polynomial<T> x4 = r.getVarPower(1, 4);
//		Polynomial<T> x6 = r.getVarPower(1, 6);
//		Polynomial<T> a = r.getEmbedding(this.a);
//		Polynomial<T> a2 = r.multiply(a, a);
//		Polynomial<T> a3 = r.multiply(a2, a);
//		Polynomial<T> b = r.getEmbedding(this.b);
//		Polynomial<T> b2 = r.multiply(b, b);
//		Polynomial<T> ab = r.multiply(a, b);
//		Polynomial<T> y2 = r.add(x3, r.multiply(a, x), b);
//		if (m == 0) {
//			this.divisionPolynomials.put(0, r.zero());
//		} else if (m == 1) {
//			this.divisionPolynomials.put(1, o);
//		} else if (m == 2) {
//			this.divisionPolynomials.put(2, o);
//		} else if (m == 3) {
//			Polynomial<T> p = r.multiply(3, x4);
//			p = r.add(p, r.multiply(6, a, x2));
//			p = r.add(p, r.multiply(12, b, x));
//			p = r.add(p, r.multiply(-1, a2));
//			this.divisionPolynomials.put(3, p);
//		} else if (m == 4) {
//			Polynomial<T> p = x6;
//			p = r.add(p, r.multiply(5, a, x4));
//			p = r.add(p, r.multiply(20, b, x3));
//			p = r.add(p, r.multiply(-5, a2, x2));
//			p = r.add(p, r.multiply(-4, ab, x));
//			p = r.add(p, r.multiply(-8, b2));
//			p = r.add(p, r.multiply(-1, a3));
//			p = r.multiply(2, p);
//			this.divisionPolynomials.put(4, p);
//		} else if (m % 2 == 0) {
//			int n = m / 2;
//			for (int i = n - 2; i <= n + 2; i++) {
//				this.getDivisionPolynomial(i);
//			}
//			Polynomial<T> psiNm2 = this.divisionPolynomials.get(n - 2);
//			Polynomial<T> psiNm1 = this.divisionPolynomials.get(n - 1);
//			Polynomial<T> psiN = this.divisionPolynomials.get(n);
//			Polynomial<T> psiN1 = this.divisionPolynomials.get(n + 1);
//			Polynomial<T> psiN2 = this.divisionPolynomials.get(n + 2);
//			Polynomial<T> p = r.subtract(r.multiply(psiN2, r.power(psiNm1, 2)), r.multiply(psiNm2, r.power(psiN1, 2)));
//			p = r.multiply(psiN, p);
//			this.divisionPolynomials.put(m, p);
//		} else if (m % 2 == 1) {
//			int n = (m - 1) / 2;
//			for (int i = n - 1; i <= n + 2; i++) {
//				this.getDivisionPolynomial(i);
//			}
//			Polynomial<T> psiNm1 = this.divisionPolynomials.get(n - 1);
//			Polynomial<T> psiN = this.divisionPolynomials.get(n);
//			Polynomial<T> psiN1 = this.divisionPolynomials.get(n + 1);
//			Polynomial<T> psiN2 = this.divisionPolynomials.get(n + 2);
//			Polynomial<T> firstTerm = r.multiply(psiN2, r.power(psiN, 3));
//			Polynomial<T> secondTerm = r.multiply(psiNm1, r.power(psiN1, 3));
//			if (n % 2 == 0) {
//				firstTerm = r.multiply(16, firstTerm, y2, y2);
//			} else {
//				secondTerm = r.multiply(16, secondTerm, y2, y2);
//			}
//			this.divisionPolynomials.put(m, r.subtract(firstTerm, secondTerm));
//		}
//		return this.getDivisionPolynomial(m);
	}

	public ProjectivePoint<T> neutral() {
		return this.pointAtInfinity;
	}

	public ProjectivePoint<T> negative(ProjectivePoint<T> p) {
		if (p.equals(pointAtInfinity)) {
			return pointAtInfinity;
		}
		return new ProjectivePoint<T>(field, p.getDehomogenisedCoord(1, 3), field.negative(
				field.add(p.getDehomogenisedCoord(2, 3), field.multiply(a1, p.getDehomogenisedCoord(1, 3)), a3)),
				field.one());
	}

	public ProjectivePoint<T> add(ProjectivePoint<T> p, ProjectivePoint<T> q) {
		return this.negative(this.getThirdIntersection(p, q));
	}

	private T computeH(ProjectivePoint<T> p, ProjectivePoint<T> q, ProjectivePoint<T> at) {
		if (p.equals(pointAtInfinity) || q.equals(pointAtInfinity)) {
			return field.one();
		}
		if (at.equals(pointAtInfinity)) {
			throw new ArithmeticException("point is in the torsion subgroup!");
		}
		T px = p.getDehomogenisedCoord(1, 3);
		T py = p.getDehomogenisedCoord(2, 3);
		T qx = q.getDehomogenisedCoord(1, 3);
		T qy = q.getDehomogenisedCoord(2, 3);
		T atx = at.getDehomogenisedCoord(1, 3);
		T aty = at.getDehomogenisedCoord(2, 3);

		T s;
		if (px.equals(qx)) {
			if (py.equals(field.zero()) || !py.equals(qy)) {
				return field.subtract(atx, px);
			}
			s = field.divide(
					field.add(field.multiply(3, px, px), field.multiply(2, a2, px), a4, field.multiply(-1, a1, py)),
					field.add(field.multiply(2, py), field.multiply(a1, px), a3));
		} else {
			s = field.divide(field.subtract(py, qy), field.subtract(px, qx));
		}
		T num = field.subtract(aty, field.add(py, field.multiply(s, field.subtract(atx, px))));
		T denom = field.subtract(field.add(atx, px, qx), field.multiply(s, s));
		return field.divide(num, denom);
	}

	private T computeF(ProjectivePoint<T> p, BigInteger n, ProjectivePoint<T> at) {
		Field<T> ff = field;
		ProjectivePoint<T> t = p;
		T f = ff.one();
		for (int i = n.bitLength() - 2; i >= 0; i--) {
			f = ff.multiply(f, f, computeH(t, t, at));
			t = add(t, t);
			if (n.testBit(i)) {
				f = ff.multiply(f, computeH(t, p, at));
				t = add(t, p);
			}
		}
		return f;
	}

	public T weilPairing(ProjectivePoint<T> p, ProjectivePoint<T> q) {
		BigInteger orderP = getOrder(p);
		BigInteger orderQ = getOrder(q);
		BigInteger gcd = orderP.gcd(orderQ);
		BigInteger n = orderP.multiply(orderQ).divide(gcd);
		return weilPairing(n, p, q);
	}

	public T weilPairing(BigInteger n, ProjectivePoint<T> p, ProjectivePoint<T> q) {
		if (!multiply(n, p).equals(neutral()) || !multiply(n, q).equals(neutral())) {
			throw new ArithmeticException("Points have wrong order!");
		}
		ProjectivePoint<T> s;
		int tries = 0;
		do {
			s = getRandomElement();
			tries++;
			if (tries > 1000) {
				throw new ArithmeticException("No points outside span!");
			}
		} while (multiply(n, s).equals(neutral()));
		ProjectivePoint<T> qs = add(q, s);
		ProjectivePoint<T> ms = negative(s);
		ProjectivePoint<T> pms = add(p, ms);
		T fpqs = computeF(p, n, qs);
		T fps = computeF(p, n, s);
		T fqpms = computeF(q, n, pms);
		T fqms = computeF(q, n, ms);
		return field.divide(field.divide(fpqs, fps), field.divide(fqpms, fqms));
	}

	@Override
	public Field<T> getField() {
		return this.field;
	}

	@SuppressWarnings("unchecked")
	@Override
	public boolean equals(Object o) {
		if (!(o instanceof EllipticCurve<?>)) {
			return false;
		}
		EllipticCurve<T> ec = (EllipticCurve<T>) o;
		return a1.equals(ec.a1) && a2.equals(ec.a2) && a3.equals(ec.a3) && a4.equals(ec.a4) && a6.equals(ec.a6);
	}

	@Override
	public int compareTo(EllipticCurve<T> o) {
		int cmp = a1.compareTo(o.a1);
		if (cmp != 0) {
			return cmp;
		}
		cmp = a2.compareTo(o.a2);
		if (cmp != 0) {
			return cmp;
		}
		cmp = a3.compareTo(o.a3);
		if (cmp != 0) {
			return cmp;
		}
		cmp = a4.compareTo(o.a4);
		if (cmp != 0) {
			return cmp;
		}
		return a6.compareTo(o.a6);
	}

	@Override
	public int getEmbeddingDimension() {
		return 2;
	}

	private T crossProductCoefficient(T x1, T y1, T x2, T y2) {
		return field.subtract(field.multiply(x1, y2), field.multiply(x2, y1));
	}

	private List<T> crossProduct(List<T> x, List<T> y) {
		List<T> result = new ArrayList<>();
		result.add(crossProductCoefficient(x.get(1), x.get(2), y.get(1), y.get(2)));
		result.add(crossProductCoefficient(x.get(2), x.get(0), y.get(2), y.get(0)));
		result.add(crossProductCoefficient(x.get(0), x.get(1), y.get(0), y.get(1)));
		return result;
	}

	public RationalFunction<T> getRationalFunction(ProjectivePoint<T> zero, ProjectivePoint<T> pole1,
			ProjectivePoint<T> pole2) {
		ProjectivePoint<T> common = this.getThirdIntersection(pole1, pole2);
		Polynomial<T> numerator;
		Polynomial<T> denominator;
		if (common.equals(zero)) {
			numerator = this.getTangentSpace(zero).get(0);
		} else {
			numerator = projectiveRing.getLinear(crossProduct(zero.getCoords(), common.getCoords()));
		}
		if (pole1.equals(pole2)) {
			denominator = this.getTangentSpace(pole1).get(0);
		} else {
			denominator = projectiveRing.getLinear(crossProduct(pole1.getCoords(), pole2.getCoords()));
		}
		return getFunctionField().getFunction(numerator, denominator);
	}

//	private List<Polynomial<T>> splitCoordinateRingElement(CoordinateRingElement<T> t) {
//		List<T> values = new ArrayList<T>();
//		values.add(null);
//		values.add(field.zero());
//		Polynomial<T> xPart = affineRing.partiallyEvaluate(t.getElement(), values);
//		values.clear();
//		values.add(null);
//		values.add(field.one());
//		Polynomial<T> yPart = affineRing.subtract(affineRing.partiallyEvaluate(t.getElement(), values), xPart);
//		List<Polynomial<T>> result = new ArrayList<>();
//		result.add(univariateRing.getEmbedding(xPart, new int[] { 0 }));
//		result.add(univariateRing.getEmbedding(yPart, new int[] { 0 }));
//		return result;
//	}
//
//	private CoordinateRingElement<T> joinCoordinateRingElement(Polynomial<T> xPart, Polynomial<T> yPart) {
//		CoordinateRing<T> cr = getCoordinateRing();
//		PolynomialRing<T> r = affineRing;
//		CoordinateRingElement<T> y = cr.getEmbedding(affineRing.getVar(2));
//		CoordinateRingElement<T> affineX = cr.getEmbedding(r.getEmbedding(xPart, new int[] { 0 }));
//		CoordinateRingElement<T> affineY = cr.getEmbedding(r.getEmbedding(yPart, new int[] { 0 }));
//		return cr.add(affineX, cr.multiply(y, affineY));
//	}

//	@Override
//	public boolean hasSimplify() {
//		return true;
//	}
//
//	@Override
//	public List<CoordinateRingElement<T>> simplify(RationalFunctionTmp<T> t) {
//		CoordinateRing<T> cr = getCoordinateRing();
//		PolynomialRing<T> r = univariateRing;
//		CoordinateRingElement<T> num = t.getNumerator();
//		CoordinateRingElement<T> denom = t.getDenominator();
//		List<Polynomial<T>> splitDenom = splitCoordinateRingElement(t.getDenominator());
//		CoordinateRingElement<T> normalizer = joinCoordinateRingElement(splitDenom.get(0),
//				r.negative(splitDenom.get(1)));
//		num = cr.multiply(normalizer, num);
//		denom = cr.multiply(normalizer, denom);
//		List<Polynomial<T>> splitNum = splitCoordinateRingElement(num);
//		Polynomial<T> numX = splitNum.get(0);
//		Polynomial<T> numY = splitNum.get(1);
//		splitDenom = splitCoordinateRingElement(denom);
//		Polynomial<T> denomX = splitDenom.get(0);
//		Polynomial<T> denomY = splitDenom.get(1);
//		if (!denomY.equals(r.zero())) {
//			throw new ArithmeticException("Could not simplify, normalizer failed!");
//		}
//		Polynomial<T> gcd = r.gcd(numX, numY);
//		gcd = r.gcd(gcd, denomX);
//		numX = r.divideChecked(numX, gcd);
//		numY = r.divideChecked(numY, gcd);
//		denomX = r.divideChecked(denomX, gcd);
//		List<CoordinateRingElement<T>> result = new ArrayList<>();
//		result.add(joinCoordinateRingElement(numX, numY));
//		result.add(joinCoordinateRingElement(denomX, r.zero()));
//		return result;
//	}

	public ProjectivePoint<T> getThirdIntersection(ProjectivePoint<T> p, ProjectivePoint<T> q) {
		T s, r1, r2;
		if (p.equals(this.pointAtInfinity))
			return this.negative(q);
		if (q.equals(this.pointAtInfinity))
			return this.negative(p);
		if (p.equals(this.negative(q)))
			return this.pointAtInfinity;
		T px = p.getDehomogenisedCoord(1, 3);
		T py = p.getDehomogenisedCoord(2, 3);
		T qx = q.getDehomogenisedCoord(1, 3);
		T qy = q.getDehomogenisedCoord(2, 3);
		if (px.equals(qx)) {
			s = field.divide(
					field.add(field.multiply(3, px, px), field.multiply(2, a2, px), field.multiply(-1, a1, py), a4),
					field.add(field.multiply(2, py), field.multiply(a1, px), a3));
		} else {
			s = this.field.divide(this.field.subtract(py, qy), this.field.subtract(px, qx));
		}
		r1 = field.multiply(s, s);
		r1 = field.subtract(r1, a2);
		r1 = field.add(r1, field.multiply(a1, s));
		r1 = field.subtract(r1, field.add(px, qx));
		r2 = field.subtract(r1, px);
		r2 = field.multiply(s, r2);
		r2 = field.add(py, r2);
//		T t = field.subtract(py, field.multiply(s, px));
//		Polynomial<T> linePolynomial = affineRing.add(affineRing.multiply(s, affineRing.getVar(1)),
//				affineRing.multiply(-1, affineRing.getVar(2)), affineRing.getEmbedding(t));
//		System.err.println(affineRing.evaluate(linePolynomial, p.getDehomogenous(3)));
//		System.err.println(affineRing.evaluate(linePolynomial, q.getDehomogenous(3)));
//		System.err.println(affineRing.evaluate(linePolynomial, new AffinePoint<T>(this.field, r1, r2)));
//		List<Polynomial<T>> intersection = new ArrayList<>();
//		intersection.add(linePolynomial);
//		intersection.add(definingPolynomial);
//		System.err.println(affineRing.solve(intersection));
//		System.err.println(p.getDehomogenous(3));
//		System.err.println(q.getDehomogenous(3));
//		System.err.println(new AffinePoint<T>(this.field, r1, r2));
//		if (!hasRationalPoint(p) || !hasRationalPoint(q)) {
//			throw new ArithmeticException("Point not on curve!");
//		}
//		if (!hasRationalPoint(new ProjectivePoint<T>(this.field, r1, r2, this.field.one()))) {
//			throw new ArithmeticException("Point not on curve!");
//		}
		return new ProjectivePoint<T>(this.field, r1, r2, this.field.one());
	}

	@Override
	public List<Polynomial<T>> getTangentSpace(ProjectivePoint<T> p) {
		List<T> line = new ArrayList<T>();
		List<Polynomial<T>> dl = getDifferentials();
		for (Polynomial<T> poly : dl)
			line.add(projectiveRing.evaluate(poly, p.getCoords()));
		return Collections.singletonList(this.projectiveRing.getLinear(line));
	}

	@Override
	public List<Polynomial<T>> getCotangentSpace(ProjectivePoint<T> p) {
		List<Polynomial<T>> list = this.getDifferentials();
		List<Polynomial<T>> reslist = new ArrayList<>();
		;
		if (!this.pointAtInfinity.equals(p)) {
			reslist.add(list.get(1));
			reslist.add(this.projectiveRing.negative(list.get(0)));
			reslist.add(this.projectiveRing.zero());
		} else {
			reslist.add(list.get(2));
			reslist.add(this.projectiveRing.zero());
			reslist.add(this.projectiveRing.negative(list.get(0)));
		}
		return reslist;
	}

	private List<Polynomial<T>> getDifferentials() {
		List<Polynomial<T>> list = new ArrayList<>();
		PolynomialRing<T> r = this.projectiveRing;
		Polynomial<T> x = r.getVar(1);
		Polynomial<T> y = r.getVar(2);
		Polynomial<T> z = r.getVar(3);
		Polynomial<T> a1 = r.getEmbedding(this.a1);
		Polynomial<T> a2 = r.getEmbedding(this.a2);
		Polynomial<T> a3 = r.getEmbedding(this.a3);
		Polynomial<T> a4 = r.getEmbedding(this.a4);
		Polynomial<T> a6 = r.getEmbedding(this.a6);
		list.add(r.add(r.multiply(3, x, x), r.multiply(2, a2, x, z), r.multiply(a4, z, z), r.multiply(-1, a1, y, z)));
		list.add(r.add(r.multiply(-2, y, z), r.multiply(-1, a1, x, z), r.multiply(-1, a3, z, z)));
		list.add(r.add(r.add(r.multiply(a2, x, x), r.multiply(2, a4, x, z), r.multiply(3, a6, z, z)),
				r.add(r.multiply(-1, y, y), r.multiply(-1, a1, x, y), r.multiply(-2, a3, y, z))));
		return list;
	}

	private List<T> getPossibleY(T x) {
		List<T> eval = new ArrayList<>();
		eval.add(x);
		eval.add(null);
		Polynomial<T> yPolynomial = affineRing.partiallyEvaluate(definingPolynomial, eval);
		UnivariatePolynomial<T> yUnivariate = univariateRing.getEmbedding(yPolynomial, new int[] { -1, 0 });
		List<T> roots = new ArrayList<>();
		roots.addAll(field.roots(yUnivariate).keySet());
		return roots;
	}

	@Override
	public ProjectivePoint<T> getRandomElement() {
		T x;
		List<T> roots;
		do {
			x = this.field.getRandomElement();
			roots = getPossibleY(x);
		} while (roots.isEmpty());
		T y = roots.get(new Random().nextInt(roots.size()));
		return new ProjectivePoint<T>(field, x, y, field.one());
	}

	private ProjectivePoint<T> getRandomElement(int degree) {
		T x;
		T y;
		int n = 0;
		BigInteger num = field.getNumberOfElements();
		BigInteger ch = field.characteristic();
		while (!num.equals(BigInteger.ONE)) {
			n++;
			num = num.divide(ch);
		}
		if (n % degree != 0) {
			throw new RuntimeException("This is wrong!");
		}
		while (true) {
			T xq = field.getRandomElement();
			x = field.zero();
			for (int i = 0; i < n; i += degree) {
				x = field.add(x, field.power(xq, ch.pow(i)));
			}
			List<T> possibleY = getPossibleY(x);
			if (possibleY.isEmpty()) {
				continue;
			}
			y = possibleY.get(new Random().nextInt(possibleY.size()));
			if (y.equals(field.power(y, ch.pow(degree)))) {
				return new ProjectivePoint<T>(field, x, y, field.one());
			}
		}
	}

	public ProjectivePoint<T> getRandomTorsionPoint(int torsion) {
		return getRandomTorsionPoint(BigInteger.valueOf(torsion));
	}

	public ProjectivePoint<T> getRandomTorsionPoint(BigInteger torsion) {
		return getRandomTorsionPoint(Integers.z().getInteger(torsion));
	}

	public ProjectivePoint<T> getRandomTorsionPoint(IntE torsion) {
		Integers z = Integers.z();
		List<ProjectivePoint<T>> basis = getTorsionPointBasis(torsion);
		ProjectivePoint<T> result = neutral();
		for (ProjectivePoint<T> basisPoint : basis) {
			result = add(result, multiply(z.getRandomElement(torsion), basisPoint));
		}
		return result;
	}

	public List<ProjectivePoint<T>> getTorsionPointBasis(int torsion) {
		return getTorsionPointBasis(BigInteger.valueOf(torsion));
	}

	public List<ProjectivePoint<T>> getTorsionPointBasis(BigInteger torsion) {
		return getTorsionPointBasis(Integers.z().getInteger(torsion));
	}

	public List<ProjectivePoint<T>> getTorsionPointBasis(IntE torsion) {
		Integers z = Integers.z();
		if (torsion.equals(z.one())) {
			return Collections.emptyList();
		}
		if (!torsionPointBasis.containsKey(torsion)) {
			IntE order = z.getInteger(getNumberOfElements());
			FactorizationResult<IntE, IntE> torsionFactors = z.uniqueFactorization(torsion);
			if (torsionFactors.primeFactors().size() > 1) {
				ProjectivePoint<T> generator1 = neutral();
				ProjectivePoint<T> generator2 = neutral();
				for (IntE prime : torsionFactors.primeFactors()) {
					IntE primePower = z.power(prime, torsionFactors.multiplicity(prime));
					List<ProjectivePoint<T>> basis = getTorsionPointBasis(primePower);
					if (basis.isEmpty()) {
						generator1 = null;
						generator2 = null;
					} else if (basis.size() == 1) {
						generator1 = generator1 != null ? add(generator1, basis.get(0)) : null;
						generator2 = null;
					} else {
						generator1 = generator1 != null ? add(generator1, basis.get(0)) : null;
						generator2 = generator2 != null ? add(generator2, basis.get(1)) : null;
					}
				}
				List<ProjectivePoint<T>> result = new ArrayList<>();
				if (generator1 != null) {
					result.add(generator1);
				}
				if (generator2 != null) {
					result.add(generator2);
				}
				torsionPointBasis.put(torsion, result);
			} else {
				if (!z.isDivisible(order, torsion)) {
					torsionPointBasis.put(torsion, Collections.emptyList());
					return Collections.emptyList();
				}
				IntE multiplier = z.divideChecked(order, torsion);
				IntE gcd;
				do {
					gcd = z.gcd(multiplier, torsion);
					multiplier = z.divideChecked(multiplier, gcd);
				} while (!gcd.equals(z.one()));
				IntE prime = torsionFactors.firstPrimeFactor();
				int multiplicity = torsionFactors.multiplicity(prime);
				IntE subPrimePower = z.power(prime, multiplicity - 1);
				int tries = 0;
				ProjectivePoint<T> torsionPoint1;
				do {
					torsionPoint1 = multiply(multiplier, getRandomElement());
					tries++;
				} while (tries < 10 && multiply(subPrimePower, torsionPoint1).equals(neutral()));
				if (multiply(subPrimePower, torsionPoint1).equals(neutral())) {
					torsionPointBasis.put(torsion, Collections.emptyList());
					return Collections.emptyList();
				}
				int multiplicationOffset = 0;
				while (!multiply(torsion, torsionPoint1).equals(neutral())) {
					multiplicationOffset++;
					torsionPoint1 = multiply(prime, torsionPoint1);
				}
				List<ProjectivePoint<T>> result = new ArrayList<>();
				result.add(torsionPoint1);
				if (z.isDivisible(order, z.power(prime, 2 * multiplicity + multiplicationOffset))) {
					ProjectivePoint<T> torsionPoint2;
					tries = 0;
					do {
						torsionPoint2 = multiply(multiplier, getRandomElement());
						tries++;
						if (multiply(subPrimePower, torsionPoint2).equals(neutral())) {
							continue;
						}
						while (!multiply(torsion, torsionPoint2).equals(neutral())) {
							torsionPoint2 = multiply(prime, torsionPoint2);
						}
					} while (tries < 10
							&& weilPairing(torsion.getValue(), torsionPoint1, torsionPoint2).equals(field.one()));
					if (tries < 10) {
						result.add(torsionPoint2);
					}
				}
				torsionPointBasis.put(torsion, result);
			}
		}
		return torsionPointBasis.get(torsion);
	}

	@Override
	public boolean isFinite() {
		return this.field.isFinite();
	}

	private void countPointsUpwards(BigInteger basecount, BigInteger q) {
		BigInteger ext = field.getNumberOfElements();
		int degree = 0;
		while (ext.mod(q).equals(BigInteger.ZERO)) {
			ext = ext.divide(q);
			degree++;
		}
		BigInteger trace = q.add(BigInteger.ONE).subtract(basecount);
		Rationals qf = Rationals.q();
		UnivariatePolynomial<Fraction> frobeniousPolynomial = qf.getUnivariatePolynomialRing()
				.getPolynomial(qf.getInteger(q), qf.getInteger(trace.negate()), qf.one());
		NumberField nf = NumberField.getNumberField(frobeniousPolynomial);
		BigInteger extendedTrace = nf.trace(nf.power(nf.alpha(), degree)).asInteger().getValue();
		this.numberOfPoints = field.getNumberOfElements().add(BigInteger.ONE).subtract(extendedTrace);
//		BigInteger a = BigInteger.TWO;
//		BigInteger b = q.add(BigInteger.ONE).subtract(basecount);
//		BigInteger base = b;
//		int n = 0;
//		BigInteger num = field.getNumberOfElements();
//		while (!num.equals(BigInteger.ONE)) {
//			n++;
//			num = num.divide(q);
//		}
//		for (int i = 1; i < n; i++) {
//			BigInteger next = base.multiply(b).subtract(q.multiply(a));
//			a = b;
//			b = next;
//		}
//		this.superSingular = true;
//		this.numberOfPoints = field.getNumberOfElements().add(BigInteger.ONE).subtract(b);
	}

	private int minimalDefinedExtension() {
		BigInteger p = field.characteristic();
		BigInteger q = field.getNumberOfElements();
		int degree = 0;
		while (!q.equals(BigInteger.ONE)) {
			degree++;
			q = q.divide(p);
		}
		q = field.getNumberOfElements();
		List<T> elements = new ArrayList<>();
		elements.add(a1);
		elements.add(a2);
		elements.add(a3);
		elements.add(a4);
		elements.add(a6);
		for (int i = 1; i <= degree; i++) {
			boolean degreeFound = true;
			for (T element : elements) {
				if (!field.power(element, p.pow(i)).equals(element)) {
					degreeFound = false;
					break;
				}
			}
			if (degreeFound) {
				return i;
			}
		}
		throw new ArithmeticException("Frobenious did not work as expected");
	}

	private boolean attemptSupersingularCount() {
		this.superSingular = false;
		this.superSingularDeterminied = true;
		if (field.characteristic().equals(BigInteger.TWO)) {
			superSingular = a1.equals(field.zero());
			return superSingular;
		}
		if (field.characteristic().equals(BigInteger.valueOf(3))) {
			superSingular = getWeierstrassForm().getRange().getA2().equals(field.zero());
			return superSingular;
		}
		T j = this.jInvariant();
		T jp = field.power(j, field.characteristic());
		T jpp = field.power(jp, field.characteristic());
		if (!j.equals(jpp)) {
			return false;
		}
		int minimalDefinedDegree = minimalDefinedExtension();
		BigInteger one = BigInteger.ONE;
		BigInteger pa1 = field.characteristic().add(one);
		BigInteger pm1 = field.characteristic().subtract(one);
		if (minimalDefinedDegree == 1) {
			if (field.characteristic().compareTo(BigInteger.valueOf(257)) <= 0) {
				int count = 1;
				for (int i = 0; i < field.characteristic().intValueExact(); i++) {
					T x = field.getInteger(i);
					for (T y : getPossibleY(x)) {
						if (field.power(y, field.characteristic()).equals(y)) {
							count++;
						}
					}
				}
				if (count == field.characteristic().intValueExact() + 1) {
					superSingular = true;
					countPointsUpwards(pa1, field.characteristic());
					return true;
				}
				return false;
			}
			// Curve defined over F_p
			for (int i = 0; i < Math.max(1, 20 - field.characteristic().bitLength() / 2); i++) {
				ProjectivePoint<T> point = this.getRandomElement(1);
				if (!multiply(pa1, point).equals(neutral())) {
					return false;
				}
			}
			superSingular = true;
			countPointsUpwards(pa1, field.characteristic());
			return true;
		}

		if (minimalDefinedDegree != 2) {
			EllipticCurve<T> newCurve = EllipticCurve.fromJInvariant(field, j);
			this.superSingular = newCurve.attemptSupersingularCount();
			this.numberOfPoints = newCurve.numberOfPoints;
			return superSingular;
		}
		// Curve defined over F_p^2
		Set<BigInteger> possibleCounts = new TreeSet<>();
		possibleCounts.add(pa1);
		possibleCounts.add(pm1);
		Set<BigInteger> remainingCounts = new TreeSet<>();
		remainingCounts.addAll(possibleCounts);
		for (int i = 0; i < Math.max(1, 20 - field.characteristic().bitLength()) || possibleCounts.size() > 1; i++) {
			for (BigInteger count : possibleCounts) {
				ProjectivePoint<T> point = getRandomElement(2);
				if (!multiply(count, point).equals(neutral())) {
					remainingCounts.remove(count);
				}
			}
			possibleCounts.clear();
			possibleCounts.addAll(remainingCounts);
		}
		if (possibleCounts.isEmpty()) {
			return false;
		}
		this.superSingular = true;
		BigInteger count = possibleCounts.iterator().next();
		countPointsUpwards(count.multiply(count), field.characteristic().pow(2));
		return true;
	}

	public void setNumberOfPointsFrom(EllipticCurve<T> domain, Isogeny<T> isogeny) {
		if (domain.equals(this) || !isogeny.getRange().equals(this) || !isogeny.getDomain().equals(domain)) {
			return;
		}
		if (domain.numberOfPoints.compareTo(BigInteger.ZERO) > 0) {
			this.numberOfPoints = domain.numberOfPoints;
		}
	}

	@Override
	public BigInteger getNumberOfElements() {
		if (!this.field.isFinite()) {
			return BigInteger.valueOf(-1);
		}
		if (this.numberOfPoints.compareTo(BigInteger.ZERO) > 0) {
			return this.numberOfPoints;
		}
		if (!isWeierstrass()) {
			this.numberOfPoints = getWeierstrassForm().getRange().getNumberOfElements();
			return this.numberOfPoints;
		}
		if (field.characteristic().equals(BigInteger.valueOf(2))) {
			return countPointsCharacteristic2();
		}
		if (field.characteristic().equals(BigInteger.valueOf(3))) {
			return countPointsCharacteristic3();
		}
		if (this.attemptSupersingularCount()) {
			return this.getNumberOfElements();
		}
		Integers z = Integers.z();
		Iterator<IntE> primeIt = z.primes();
		BigInteger zero = BigInteger.ZERO;
		BigInteger one = BigInteger.ONE;
		BigInteger two = BigInteger.TWO;
		List<BigInteger> primes = new ArrayList<BigInteger>();
		List<BigInteger> t = new ArrayList<BigInteger>();
		IntE i = primeIt.next();
		BigInteger product = one;
		BigInteger q = this.field.getNumberOfElements();
		BigInteger twoSqrtQ = q.sqrt().add(one).shiftLeft(1);
		while (product.compareTo(twoSqrtQ.shiftLeft(1)) < 0) {
			if (i.getValue().mod(this.field.characteristic()).equals(zero)) {
				i = primeIt.next();
				continue;
			}
			if (i.getValue().subtract(one).mod(this.field.characteristic()).equals(zero)) {
				i = primeIt.next();
				continue;
			}
			if (i.getValue().add(one).mod(this.field.characteristic()).equals(zero)) {
				i = primeIt.next();
				continue;
			}
			primes.add(i.getValue());
			product = product.multiply(i.getValue());
			i = primeIt.next();
		}
		PolynomialRing<T> r = this.affineRing;
		Polynomial<T> x = r.getVar(1);
		Polynomial<T> y = r.getVar(2);
		for (BigInteger l : primes) {
			if (l.equals(two)) {
				// Computing whether curve has two torsion points.
				if (field.hasRoots(getAdjustedRhs())) {
					t.add(zero);
				} else {
					t.add(one);
				}
				continue;
			}
			int qbar = q.mod(l).intValue();
			if ((l.intValue() - 1) / 2 > 0) {
				qbar -= l.intValue();
			}
			List<Polynomial<T>> idealGen = new ArrayList<Polynomial<T>>();
			Polynomial<T> psiL = this.getDivisionPolynomial(l.intValue());
			idealGen.add(r.getEmbedding(psiL, new int[] { 0, 1 }));
			idealGen.add(r.getEmbedding(this.definingPolynomial, new int[] { 0, 1 }));
			PolynomialIdeal<T> ideal = r.getIdeal(idealGen);
			psiL = this.univariateRing.getEmbedding(psiL, new int[] { 0 });
			CoordinateRing<T> cr = ideal.divideOut();
			// Calculating x^q, y^q, x^(q^2), y^(q^2)
			CoordinateRingElement<T> xq = cr.power(cr.getEmbedding(x), q);
			CoordinateRingElement<T> yq = cr.power(cr.getEmbedding(y), q);
			List<CoordinateRingElement<T>> xyq = new ArrayList<>();
			xyq.add(xq);
			xyq.add(yq);
			CoordinateRingElement<T> xqq = cr.power(xq, q);
			CoordinateRingElement<T> yqq = cr.power(yq, q);
			List<CoordinateRingElement<T>> qxy = multiplyGeneric(qbar, psiL, cr);
			CoordinateRingElement<T> qx = qxy.get(0);
			CoordinateRingElement<T> qy = qxy.get(1);
			if (!xqq.equals(qx)) {
				List<CoordinateRingElement<T>> lhs = this.addGeneric(xqq, yqq, qx, qy, psiL, cr);
				CoordinateRingElement<T> lhsX = lhs.get(0);
				CoordinateRingElement<T> lhsY = lhs.get(1);
				boolean found = false;
				for (int tc = 1; tc <= (l.intValue() - 1) / 2; tc++) {
					List<CoordinateRingElement<T>> candidate = this.multiplyGeneric(tc, psiL, cr);
					CoordinateRingElement<T> tcqx = cr.power(candidate.get(0), q);// cr.substitute(candidate.get(0),
																					// xyq);
					if (tcqx.equals(lhsX)) {
						found = true;
						CoordinateRingElement<T> tcqy = cr.power(candidate.get(1), q);// cr.substitute(candidate.get(1),
																						// xyq);
						if (tcqy.equals(lhsY)) {
							t.add(BigInteger.valueOf(tc));
							break;
						} else if (tcqy.equals(cr.negative(lhsY))) {
							t.add(BigInteger.valueOf(-tc));
							break;
						} else {
							throw new ArithmeticException("No y coordinate found.");
						}
					}
				}
				if (!found) {
					throw new ArithmeticException("No x coordinate found.");
				}
			} else {
				PrimeField lField = PrimeField.getPrimeField(l);
				PFE qModL = lField.getInteger(q);
				if (lField.hasSqrt(qModL)) {
					PFE w = lField.sqrt(qModL).keySet().iterator().next();
					List<CoordinateRingElement<T>> candidate = this.multiplyGeneric(w.getValue().intValueExact(), psiL,
							cr);
					CoordinateRingElement<T> wqx = cr.substitute(candidate.get(0), xyq);
					if (wqx.equals(xqq)) {
						CoordinateRingElement<T> wqy = cr.substitute(candidate.get(1), xyq);
						if (wqy.equals(yqq)) {
							t.add(w.getValue().shiftLeft(1));
						} else if (wqy.equals(cr.negative(yqq))) {
							t.add(w.getValue().shiftLeft(1).negate());
						} else {
							t.add(zero);
						}
					} else {
						t.add(zero);
					}
				} else {
					t.add(zero);
				}
			}
		}
		BigInteger result = MiscAlgorithms.chineseRemainder(t, primes, true);
		this.numberOfPoints = this.field.getNumberOfElements().add(one).subtract(result);
		return this.getNumberOfElements();
	}

	private Polynomial<T> invertXOnlyPolynomial(Polynomial<T> a, Polynomial<T> psiL) {
		PolynomialRing<T> univarring = this.univariateRing;
		Polynomial<T> univarA = univarring.getEmbedding(a, new int[] { 0 });
		// Inverting psiQsq mod psiL:
		ExtendedEuclideanResult<Polynomial<T>> egcd = univarring.extendedEuclidean(univarA, psiL);
		Polynomial<T> univarAInv = univarring.multiply(egcd.getCoeff1(), univarring.inverse(egcd.getGcd()));
		return this.affineRing.getEmbedding(univarAInv, new int[] { 0 });
	}

	private List<CoordinateRingElement<T>> addGeneric(CoordinateRingElement<T> x1, CoordinateRingElement<T> y1,
			CoordinateRingElement<T> x2, CoordinateRingElement<T> y2, Polynomial<T> psiL, CoordinateRing<T> cr) {
		PolynomialRing<T> r = this.univariateRing;
		CoordinateRingElement<T> ydiff = cr.subtract(y2, y1);
		CoordinateRingElement<T> xdiff = cr.subtract(x2, x1);
		Polynomial<T> xdiffPoly = r.getEmbedding(xdiff.getElement(), new int[] { 0 });
		Polynomial<T> gcdPoly = r.gcd(xdiffPoly, psiL);
		Polynomial<T> gcdMulti = this.affineRing.getEmbedding(gcdPoly, new int[] { 0 });
		xdiffPoly = r.divideChecked(xdiffPoly, gcdPoly);
		CoordinateRingElement<T> ydiffReduced = cr
				.getEmbedding(this.affineRing.divideChecked(ydiff.getElement(), gcdMulti));
		CoordinateRingElement<T> xdiffInvReduced = cr.getEmbedding(this.invertXOnlyPolynomial(xdiffPoly, psiL));
		CoordinateRingElement<T> s = cr.multiply(ydiffReduced, xdiffInvReduced);
		CoordinateRingElement<T> xr = cr.add(cr.multiply(s, s), cr.multiply(a1, s),
				cr.getEmbedding(field.multiply(-1, a2)));
		xr = cr.subtract(xr, cr.add(x1, x2));
		CoordinateRingElement<T> yr = cr.add(y1, cr.multiply(s, cr.subtract(xr, x1)));
		List<CoordinateRingElement<T>> result = new ArrayList<>();
		result.add(xr);
		result.add(cr.add(cr.negative(yr), cr.negative(cr.multiply(a1, xr)), cr.negative(cr.getEmbedding(a3))));
		return result;
	}

	private List<CoordinateRingElement<T>> multiplyGeneric(int n, Polynomial<T> psiL, CoordinateRing<T> cr) {
		// Calculating [qbar] (x, y)
		PolynomialRing<T> r = this.affineRing;
		int nabs = Math.abs(n);
		Polynomial<T> psiNsq = this.getDivisionPolynomial(nabs);
		psiNsq = r.multiply(psiNsq, psiNsq);
		psiNsq = cr.getPolynomialRing()
				.getEmbedding(r.reduce(psiNsq, Collections.singletonList(this.definingPolynomial)), new int[] { 0, 1 });
		CoordinateRingElement<T> psiNsqInv = cr.getEmbedding(this.invertXOnlyPolynomial(psiNsq, psiL));
		CoordinateRingElement<T> psiNm1 = cr.getEmbedding(this.getDivisionPolynomial(nabs == 0 ? 1 : (nabs - 1)));
		CoordinateRingElement<T> psiNp1 = cr.getEmbedding(this.getDivisionPolynomial(nabs + 1));
		CoordinateRingElement<T> qx = cr.subtract(cr.getEmbedding(r.getVar(1)), cr.multiply(psiNsqInv, psiNm1, psiNp1));
		CoordinateRingElement<T> qy = cr.multiply(cr.getEmbedding(this.getDivisionPolynomial(2 * nabs)), psiNsqInv,
				psiNsqInv);
		qy = cr.multiply(qy, cr.getEmbedding(r.getEmbedding(this.field.inverse(this.field.getInteger(2)))));
		if (n < 0) {
			qy = cr.negative(qy);
		}
		List<CoordinateRingElement<T>> result = new ArrayList<>();
		result.add(qx);
		result.add(qy);
		return result;
	}

	private Polynomial<T> secondModularPolynomial() {
		Integers z = Integers.z();
		PolynomialRing<T> r = affineRing;
		Polynomial<T> phi2 = r.add(r.getVarPower(1, 3), r.getVarPower(1, 3));
		phi2 = r.subtract(phi2, r.multiply(r.getVarPower(1, 2), r.getVarPower(2, 2)));
		phi2 = r.add(phi2, r.multiply(2 * 2 * 2 * 2 * 3 * 31,
				r.add(r.multiply(r.getVarPower(1, 2), r.getVar(2)), r.multiply(r.getVar(1), r.getVarPower(2, 2)))));
		phi2 = r.subtract(phi2, r.multiply(
				z.multiply(z.power(z.getInteger(2), 4), z.power(z.getInteger(3), 4), z.power(z.getInteger(5), 3)),
				r.add(r.getVarPower(1, 2), r.getVarPower(2, 2))));
		phi2 = r.add(phi2,
				r.multiply(z.multiply(z.power(z.getInteger(3), 4), z.power(z.getInteger(5), 3), z.getInteger(4027)),
						r.getVar(1), r.getVar(2)));
		phi2 = r.add(phi2, r.multiply(
				z.multiply(z.power(z.getInteger(2), 8), z.power(z.getInteger(3), 7), z.power(z.getInteger(5), 6)),
				r.add(r.getVar(1), r.getVar(2))));
		phi2 = r.subtract(phi2, r.getInteger(
				z.multiply(z.power(z.getInteger(2), 12), z.power(z.getInteger(3), 9), z.power(z.getInteger(5), 9))));
		return phi2;
	}

	private BigInteger countPointsCharacteristic2() {
		this.superSingularDeterminied = true;
		if (a1.equals(field.zero())) {
			this.superSingular = true;
		}
		int degree = minimalDefinedExtension();
		if (degree == 1) {
			int count = 1;
			for (int i = 0; i < 2; i++) {
				T x = field.getInteger(i);
				List<T> ys = getPossibleY(x);
				for (T y : ys) {
					if (field.power(y, 2).equals(y)) {
						count++;
					}
				}
			}
			countPointsUpwards(BigInteger.valueOf(count), BigInteger.TWO);
			return numberOfPoints;
		} else if (degree == 2) {
			FiniteField f4 = FiniteField.getFiniteField(4);
			FieldEmbedding<PFE, FFE, FiniteField> embedding = new FieldEmbedding<>((FiniteField) field, f4);
			int count = 1;
			for (FFE xBase : f4) {
				@SuppressWarnings("unchecked")
				T x = (T) embedding.getEmbedding(xBase);
				List<T> ys = getPossibleY(x);
				for (T y : ys) {
					if (field.power(y, 4).equals(y)) {
						count++;
					}
				}
			}
			countPointsUpwards(BigInteger.valueOf(count), BigInteger.valueOf(4));
			return numberOfPoints;
		}
		if (superSingular) {
			Rationals q = Rationals.q();
			UnivariatePolynomialRing<Fraction> rationalPolynomials = q.getUnivariatePolynomialRing();
			Fraction one = q.one();
			Fraction two = q.getInteger(2);
			Fraction four = q.getInteger(2);
			Set<IntE> possibleTraces = new TreeSet<>();
			NumberField nf = NumberField.getNumberField(rationalPolynomials.getPolynomial(two, two, one));
			possibleTraces.add(nf.trace(nf.power(nf.alpha(), degree)).asInteger());
			nf = NumberField.getNumberField(rationalPolynomials.getPolynomial(two, q.zero(), one));
			possibleTraces.add(nf.trace(nf.power(nf.alpha(), degree)).asInteger());
			nf = NumberField.getNumberField(rationalPolynomials.getPolynomial(two, q.negative(two), one));
			possibleTraces.add(nf.trace(nf.power(nf.alpha(), degree)).asInteger());
			if (degree % 2 == 0) {
				nf = NumberField.getNumberField(rationalPolynomials.getPolynomial(four, two, one));
				possibleTraces.add(nf.trace(nf.power(nf.alpha(), degree / 2)).asInteger());
				nf = NumberField.getNumberField(rationalPolynomials.getPolynomial(four, q.zero(), one));
				possibleTraces.add(nf.trace(nf.power(nf.alpha(), degree / 2)).asInteger());
				nf = NumberField.getNumberField(rationalPolynomials.getPolynomial(four, q.negative(two), one));
				possibleTraces.add(nf.trace(nf.power(nf.alpha(), degree / 2)).asInteger());
			}
			BigInteger p1 = field.characteristic().pow(degree).add(BigInteger.ONE);
			Set<IntE> remainingTraces = new TreeSet<>();
			remainingTraces.addAll(possibleTraces);
			while (possibleTraces.size() > 1) {
				for (IntE trace : possibleTraces) {
					ProjectivePoint<T> point = getRandomElement(degree);
					BigInteger number = p1.subtract(trace.getValue());
					if (!multiply(number, point).equals(neutral())) {
						remainingTraces.remove(trace);
					}
				}
				possibleTraces.clear();
				possibleTraces.addAll(remainingTraces);
			}
			countPointsUpwards(p1.subtract(possibleTraces.iterator().next().getValue()),
					field.characteristic().pow(degree));
			return numberOfPoints;
		}
		T traceA2 = finiteFieldTrace(a2);
		if (traceA2.equals(field.zero()) && !a2.equals(field.zero())) {
			List<T> coeffs = new ArrayList<>();
			coeffs.add(a2);
			coeffs.add(field.one());
			coeffs.add(field.one());
			T xyOffset = field.roots(univariateRing.getPolynomial(coeffs)).keySet().iterator().next();
			EllipticCurve<T> reduced = new Isomorphism<>(this, field.one(), field.zero(), xyOffset, field.zero())
					.getRange();
			numberOfPoints = reduced.getNumberOfElements();
			return numberOfPoints;
		} else if (traceA2.equals(field.one())) {
			numberOfPoints = field.getNumberOfElements().multiply(BigInteger.TWO).add(BigInteger.TWO)
					.subtract(getQuadraticTwist().getNumberOfElements());
			return numberOfPoints;
		}
		return null;
	}

	private BigInteger countPointsCharacteristic3() {
		this.superSingularDeterminied = true;
		if (a2.equals(field.zero())) {
			this.superSingular = true;
		}
		int degree = minimalDefinedExtension();

		if (superSingular) {
			Rationals q = Rationals.q();
			UnivariatePolynomialRing<Fraction> rationalPolynomials = q.getUnivariatePolynomialRing();
			Fraction one = q.one();
			Fraction two = q.getInteger(2);
			Fraction four = q.getInteger(2);
			Set<IntE> possibleTraces = new TreeSet<>();
			NumberField nf = NumberField.getNumberField(rationalPolynomials.getPolynomial(two, two, one));
			possibleTraces.add(nf.trace(nf.power(nf.alpha(), degree)).asInteger());
			nf = NumberField.getNumberField(rationalPolynomials.getPolynomial(two, q.zero(), one));
			possibleTraces.add(nf.trace(nf.power(nf.alpha(), degree)).asInteger());
			nf = NumberField.getNumberField(rationalPolynomials.getPolynomial(two, q.negative(two), one));
			possibleTraces.add(nf.trace(nf.power(nf.alpha(), degree)).asInteger());
			if (degree % 2 == 0) {
				nf = NumberField.getNumberField(rationalPolynomials.getPolynomial(four, two, one));
				possibleTraces.add(nf.trace(nf.power(nf.alpha(), degree / 2)).asInteger());
				nf = NumberField.getNumberField(rationalPolynomials.getPolynomial(four, q.zero(), one));
				possibleTraces.add(nf.trace(nf.power(nf.alpha(), degree / 2)).asInteger());
				nf = NumberField.getNumberField(rationalPolynomials.getPolynomial(four, q.negative(two), one));
				possibleTraces.add(nf.trace(nf.power(nf.alpha(), degree / 2)).asInteger());
			}
			BigInteger p1 = field.characteristic().pow(degree).add(BigInteger.ONE);
			Set<IntE> remainingTraces = new TreeSet<>();
			remainingTraces.addAll(possibleTraces);
			while (possibleTraces.size() > 1) {
				for (IntE trace : possibleTraces) {
					ProjectivePoint<T> point = getRandomElement(degree);
					BigInteger number = p1.subtract(trace.getValue());
					if (!multiply(number, point).equals(neutral())) {
						remainingTraces.remove(trace);
					}
				}
				possibleTraces.clear();
				possibleTraces.addAll(remainingTraces);
			}
			countPointsUpwards(p1.subtract(possibleTraces.iterator().next().getValue()),
					field.characteristic().pow(degree));
			return numberOfPoints;
		}
		return null;
	}

	public ProjectiveMorphism<T> translationMorphism(ProjectivePoint<T> point) {
		return null;
//		if (point.equals(pointAtInfinity)) {
//			List<Polynomial<T>> values = new ArrayList<>();
//			values.add(projectiveRing.getVar(1));
//			values.add(projectiveRing.getVar(2));
//			values.add(projectiveRing.getVar(3));
//			return new ProjectiveMorphism<>(asGenericProjectiveScheme(), asGenericProjectiveScheme(), values);
//		}
//		PolynomialRing<T> r = projectiveRing;
//		Polynomial<T> xp = r.getEmbedding(point.getDehomogenisedCoord(1, 3));
//		Polynomial<T> yp = r.getEmbedding(point.getDehomogenisedCoord(2, 3));
//		Polynomial<T> xmxp = r.subtract(r.getVar(1), r.multiply(xp, r.getVar(3)));
//		Polynomial<T> ymyp = r.subtract(r.getVar(2), r.multiply(yp, r.getVar(3)));
//		Polynomial<T> x2 = r.getVarPower(1, 2);
//		Polynomial<T> xz = r.multiply(r.getVar(1), r.getVar(3));
//		Polynomial<T> yz = r.multiply(r.getVar(2), r.getVar(3));
//		Polynomial<T> z2 = r.getVarPower(3, 2);
//		Polynomial<T> xpx2 = r.multiply(xp, x2);
//		Polynomial<T> xp2axz = r.multiply(r.add(r.power(xp, 2), r.getEmbedding(a)), xz);
//		Polynomial<T> yyz = r.multiply(-2, yp, yz);
//		Polynomial<T> cz2 = r.multiply(r.add(r.multiply(2, r.getEmbedding(b)), r.multiply(a, xp)), z2);
//		Polynomial<T> xR = r.add(xpx2, xp2axz, yyz, cz2);
//		Polynomial<T> xResult = r.multiply(xmxp, xR);
//		Polynomial<T> yResult = r.negative(r.add(r.multiply(yp, r.power(xmxp, 3)),
//				r.multiply(ymyp, r.subtract(xR, r.multiply(xp, r.power(xmxp, 2))))));
//		Polynomial<T> zResult = r.power(xmxp, 3);
//		List<Polynomial<T>> values = new ArrayList<>();
//		values.add(xResult);
//		values.add(yResult);
//		values.add(zResult);
//		return new ProjectiveMorphism<>(asGenericProjectiveScheme(), asGenericProjectiveScheme(), values);
	}

	public ProjectiveMorphism<T> multiplicationMorphism(int n) {
		boolean negative = false;
		if (n < 0) {
			negative = true;
			n = -n;
		}
		if (n == 0) {
			List<Polynomial<T>> asPolynomials = new ArrayList<>();
			asPolynomials.add(projectiveRing.getEmbedding(field.zero()));
			asPolynomials.add(projectiveRing.getEmbedding(field.one()));
			asPolynomials.add(projectiveRing.getEmbedding(field.zero()));
			return new ProjectiveMorphism<>(asGenericProjectiveScheme(), asGenericProjectiveScheme(), asPolynomials);
		}
		if (n == 1) {
			return new Isomorphism<>(this, negative ? field.getInteger(-1) : field.getInteger(1), field.zero(),
					field.zero(), field.zero()).asMorphism();
		}
		CoordinateRing<T> r = getCoordinateRing();
		CoordinateRingElement<T> psiNm1 = r.getEmbedding(getDivisionPolynomial(n - 1));
		CoordinateRingElement<T> psiN = r.getEmbedding(getDivisionPolynomial(n));
		CoordinateRingElement<T> psiNp1 = r.getEmbedding(getDivisionPolynomial(n + 1));
		Polynomial<T> phiN = r
				.subtract(r.multiply(r.getVar(1), r.power(psiN, 2)), coordinateRing.multiply(psiNm1, psiNp1))
				.getElement();
		Polynomial<T> omegaN = affineRing.divideChecked(affineRing.subtract(
				affineRing.multiply(getDivisionPolynomial(n + 2), getDivisionPolynomial(n - 1),
						getDivisionPolynomial(n - 1)),
				affineRing.multiply(getDivisionPolynomial(n - 2), getDivisionPolynomial(n + 1),
						getDivisionPolynomial(n + 1))),
				affineRing.multiply(4, affineRing.getVar(2)));
		List<Polynomial<T>> polynomials = new ArrayList<>();
		Polynomial<T> psiN1 = psiN.getElement();
		Polynomial<T> psiN3 = r.multiply(psiN, psiN, psiN).getElement();
		Polynomial<T> x = projectiveRing.multiply(affineRing.homogenize(phiN), affineRing.homogenize(psiN1));
		Polynomial<T> y = affineRing.homogenize(!negative ? omegaN : affineRing.negative(omegaN));
		Polynomial<T> z = affineRing.homogenize(psiN3);
		int degreeMax = Math.max(Math.max(x.degree(), y.degree()), z.degree());
		x = projectiveRing.multiply(x, projectiveRing.power(projectiveRing.getVar(3), degreeMax - x.degree()));
		y = projectiveRing.multiply(y, projectiveRing.power(projectiveRing.getVar(3), degreeMax - y.degree()));
		z = projectiveRing.multiply(z, projectiveRing.power(projectiveRing.getVar(3), degreeMax - z.degree()));
		polynomials.add(x);
		polynomials.add(y);
		polynomials.add(z);
		return new ProjectiveMorphism<>(asGenericProjectiveScheme(), asGenericProjectiveScheme(), polynomials);
	}

	private UnivariatePolynomial<T> getAdjustedRhs() {
		if (field.characteristic().equals(BigInteger.valueOf(2))) {
			throw new ArithmeticException("Not defined!");
		}
		if (!isWeierstrass()) {
			return getWeierstrassForm().getRange().getAdjustedRhs();
		}
		List<T> coefficients = new ArrayList<>();
		coefficients.add(a6);
		coefficients.add(a4);
		coefficients.add(a2);
		coefficients.add(field.one());
		return univariateRing.getPolynomial(coefficients);
	}

	public List<ProjectivePoint<T>> getTorsionPoints(int l) {
		if (l == 1) {
			return Collections.singletonList(this.pointAtInfinity);
		}
		List<ProjectivePoint<T>> torsionPoints = new ArrayList<>(l * l);
		if (l == 2) {
			torsionPoints.add(this.pointAtInfinity);
			if (field.characteristic().equals(BigInteger.TWO)) {
				if (!a1.equals(field.zero())) {
					T x = field.divide(a3, a1);
					T y = getPossibleY(x).get(0);
					torsionPoints.add(new ProjectivePoint<>(field, x, y, field.one()));
				}
			} else {
				Map<T, Integer> rhsRoots = this.field.roots(getAdjustedRhs());
				for (T root : rhsRoots.keySet()) {
					torsionPoints
							.add(new ProjectivePoint<T>(this.field, root, getPossibleY(root).get(0), this.field.one()));
				}
			}
		}
		if (l > 2) {
			Polynomial<T> divPoly = null;
			if (l % 2 == 0) {
				torsionPoints.addAll(this.getTorsionPoints(2));
				divPoly = this.univariateRing.getEmbedding(
						univariateRing.divideChecked(getDivisionPolynomial(l), getDivisionPolynomial(2)),
						new int[] { 0 });
			} else {
				this.getDivisionPolynomial(l);
				divPoly = this.univariateRing.getEmbedding(this.divisionPolynomials.get(l), new int[] { 0 });
				torsionPoints.add(this.pointAtInfinity);
			}
			Map<T, Integer> xs = this.field.roots(divPoly);
			for (T x : xs.keySet()) {
				for (T y : this.field.sqrt(this.univariateRing.evaluate(this.rhs, Collections.singletonList(x)))
						.keySet()) {
					torsionPoints.add(new ProjectivePoint<T>(this.field, x, y, this.field.one()));
				}
			}
		}
		return torsionPoints;
	}

	@Override
	public Iterator<ProjectivePoint<T>> iterator() {
		return new Iterator<ProjectivePoint<T>>() {
			private ProjectivePoint<T> nextPoint = pointAtInfinity;
			private ProjectivePoint<T> nextNextPoint = null;
			private Iterator<T> xIterator = field.iterator();

			@Override
			public boolean hasNext() {
				return nextPoint != null;
			}

			@Override
			public ProjectivePoint<T> next() {
				ProjectivePoint<T> point = nextPoint;
				if (nextNextPoint != null) {
					nextPoint = nextNextPoint;
					nextNextPoint = null;
					return point;
				}
				T x;
				List<T> ys;
				do {
					if (!xIterator.hasNext()) {
						nextPoint = null;
						return point;
					}
					x = xIterator.next();
					ys = getPossibleY(x);
				} while (ys.isEmpty());
				Iterator<T> yIterator = ys.iterator();
				nextPoint = new ProjectivePoint<T>(field, x, yIterator.next(), field.one());
				if (yIterator.hasNext()) {
					nextNextPoint = new ProjectivePoint<T>(field, x, yIterator.next(), field.one());
				}
				return point;
			}
		};
	}

	@Override
	public ProjectivePoint<T> inverse(ProjectivePoint<T> t) {
		return this.negative(t);
	}

	@Override
	public ProjectivePoint<T> operate(ProjectivePoint<T> t1, ProjectivePoint<T> t2) {
		return this.add(t1, t2);
	}

	public ProjectivePoint<T> multiply(int n, ProjectivePoint<T> t) {
		return this.multiply(BigInteger.valueOf(n), t);
	}

	public ProjectivePoint<T> multiply(BigInteger n, ProjectivePoint<T> t) {
		return this.power(n, t);
	}

	public ProjectivePoint<T> multiply(IntE n, ProjectivePoint<T> t) {
		return this.multiply(n.getValue(), t);
	}

	@Override
	public String toString() {
		return lhs + " = " + rhs;
	}

	@Override
	public List<RationalFunction<T>> getRiemannRochSpace(WeilDivisor<T> div) {
		if (div.getDegree() < 0)
			return Collections.emptyList();
		List<RationalFunction<T>> functions = new ArrayList<>();
		FunctionField<T> ff = getFunctionField();
		List<ProjectivePoint<T>> zeroes = div.getPoles();
		List<ProjectivePoint<T>> poles = div.getZeroes();
		if (poles.size() == 0) {
			return Collections.singletonList(ff.one());
		}
		RationalFunction<T> f = ff.one();
		ProjectivePoint<T> firstpole = poles.get(0);
		ProjectivePoint<T> lastzero = null;
		int size = zeroes.size();
		if (zeroes.size() == poles.size()) {
			size--;
			lastzero = zeroes.get(size);
		}
		for (int i = 0; i < size; i++) {
			ProjectivePoint<T> zero = zeroes.get(i);
			ProjectivePoint<T> pole = poles.get(i + 1);
			f = ff.multiply(f, this.getRationalFunction(zero, firstpole, pole));
			firstpole = this.getThirdIntersection(firstpole, pole);
			firstpole = this.getThirdIntersection(firstpole, zero);
		}
		if (lastzero != null) {
			if (firstpole.equals(lastzero))
				return Collections.singletonList(f);
			else
				return Collections.emptyList();
		}
		functions.add(f);
		for (int i = size + 1; i < poles.size(); i++) {
			ProjectivePoint<T> pole = poles.get(i);
			ProjectivePoint<T> zero = this.getThirdIntersection(firstpole, pole);
			while (pole.equals(zero) || firstpole.equals(zero))
				zero = this.getRandomElement();
			functions.add(ff.multiply(f, this.getRationalFunction(zero, firstpole, pole)));
		}
		return functions;
	}

	@Override
	public boolean isPrincipal(WeilDivisor<T> div) {
		if (div.getDegree() != 0)
			return false;
		return this.getRiemannRochSpace(div).size() == 1;
	}

	@Override
	public List<EllipticCurve<T>> irreducibleComponents() {
		return Collections.singletonList(this);
	}

	@Override
	public EllipticCurve<T> reduced() {
		return this;
	}
}
