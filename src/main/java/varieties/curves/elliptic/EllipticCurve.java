package varieties.curves.elliptic;

import java.lang.reflect.Array;
import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Deque;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.Random;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;

import fields.finitefields.FiniteField;
import fields.finitefields.FiniteField.FFE;
import fields.finitefields.ModuloIntegerRing;
import fields.finitefields.ModuloIntegerRing.ModuloIntegerRingElement;
import fields.finitefields.PrimeField;
import fields.finitefields.PrimeField.PFE;
import fields.floatingpoint.Complex;
import fields.floatingpoint.Complex.ComplexNumber;
import fields.floatingpoint.Reals;
import fields.floatingpoint.Reals.Real;
import fields.helper.AbstractElement;
import fields.helper.FieldEmbedding;
import fields.integers.Integers;
import fields.integers.Integers.IntE;
import fields.integers.Rationals;
import fields.integers.Rationals.Fraction;
import fields.interfaces.AlgebraicExtensionElement;
import fields.interfaces.DedekindRing;
import fields.interfaces.DedekindRingExtension;
import fields.interfaces.DiscreteValuationField;
import fields.interfaces.DiscreteValuationFieldExtension;
import fields.interfaces.Element;
import fields.interfaces.Field;
import fields.interfaces.Field.Extension;
import fields.interfaces.FieldExtension;
import fields.interfaces.FieldExtension.SplittingFieldResult;
import fields.interfaces.GlobalField;
import fields.interfaces.GlobalField.ExtensionOfGlobalField;
import fields.interfaces.GlobalFieldExtension;
import fields.interfaces.Group;
import fields.interfaces.Ideal;
import fields.interfaces.MathMap;
import fields.interfaces.Polynomial;
import fields.interfaces.PolynomialRing;
import fields.interfaces.Ring.FactorizationResult;
import fields.interfaces.Ring.QuotientAndRemainderResult;
import fields.interfaces.UnivariatePolynomial;
import fields.interfaces.UnivariatePolynomialRing;
import fields.interfaces.ValueField;
import fields.local.PAdicField;
import fields.local.PAdicField.PAdicNumber;
import fields.local.Value;
import fields.numberfields.CompletedNumberField;
import fields.numberfields.CompletedNumberField.Ext;
import fields.numberfields.EmbeddedNumberField;
import fields.numberfields.LocalizedNumberField;
import fields.numberfields.NumberField;
import fields.numberfields.NumberField.NFE;
import fields.numberfields.NumberFieldIntegers;
import fields.numberfields.NumberFieldIntegers.NumberFieldIdeal;
import fields.numberfields.NumberFieldOrder;
import fields.numberfields.PicardGroup.OrderIdealClass;
import fields.polynomials.AbstractPolynomialRing;
import fields.polynomials.CoordinateRing;
import fields.polynomials.CoordinateRing.CoordinateRingElement;
import fields.polynomials.Monomial;
import fields.polynomials.PolynomialIdeal;
import fields.vectors.FreeModule;
import fields.vectors.FreeSubModule;
import fields.vectors.Matrix;
import fields.vectors.Vector;
import util.FunctionMathMap;
import util.MiscAlgorithms;
import util.Pair;
import util.graph.Graph;
import util.graph.Graph.GraphBuilder;
import varieties.FunctionField;
import varieties.RationalFunction;
import varieties.affine.AffinePoint;
import varieties.curves.ProjectiveLine;
import varieties.curves.SmoothCurve;
import varieties.curves.WeilDivisor;
import varieties.projective.AbstractProjectiveScheme;
import varieties.projective.GenericProjectiveScheme;
import varieties.projective.ProjectiveMorphism;
import varieties.projective.ProjectivePoint;

public class EllipticCurve<T extends Element<T>> extends AbstractProjectiveScheme<T>
		implements SmoothCurve<T>, Group<ProjectivePoint<T>>, Element<EllipticCurve<T>> {

	public static class Isomorphism<T extends Element<T>> extends AbstractElement<Isogeny<T>> implements Isogeny<T> {
		private final EllipticCurve<T> domain;
		private final EllipticCurve<T> range;
		private final Field<T> field;
		private final T scale;
		private final T xOffset;
		private final T xyOffset;
		private final T yOffset;

		public Isomorphism(EllipticCurve<T> domain, T scale, T xOffset, T xyOffset, T yOffset) {
			this.domain = domain;
			this.field = domain.getField();
			this.scale = scale;
			this.xOffset = xOffset;
			this.xyOffset = xyOffset;
			this.yOffset = yOffset;
			T a1 = domain.getA1();
			T a2 = domain.getA2();
			T a3 = domain.getA3();
			T a4 = domain.getA4();
			T a6 = domain.getA6();
			T a1Range = field.divide(field.add(a1, field.multiply(2, xyOffset)), scale);
			T a2Range = field.divide(field.add(a2, field.multiply(-1, xyOffset, a1), field.multiply(3, xOffset),
					field.multiply(-1, xyOffset, xyOffset)), field.power(scale, 2));
			T a3Range = field.divide(field.add(a3, field.multiply(xOffset, a1), field.multiply(2, yOffset)),
					field.power(scale, 3));
			T a4Range = field.divide(
					field.add(a4, field.add(field.multiply(-1, xyOffset, a3), field.multiply(2, xOffset, a2)),
							field.add(field.multiply(-1, field.add(yOffset, field.multiply(xOffset, xyOffset)), a1),
									field.multiply(3, xOffset, xOffset), field.multiply(-2, xyOffset, yOffset))),
					field.power(scale, 4));
			T a6Range = field.divide(field.add(a6,
					field.add(field.multiply(xOffset, a4), field.multiply(xOffset, xOffset, a2),
							field.power(xOffset, 3)),
					field.add(field.multiply(-1, yOffset, a3), field.multiply(-1, yOffset, yOffset),
							field.multiply(-1, xOffset, yOffset, a1))),
					field.power(scale, 6));
			this.range = new EllipticCurve<>(field, a1Range, a2Range, a3Range, a4Range, a6Range);// fromScaledPolynomial(substituted);
		}

		@Override
		public String toString() {
			Isomorphism<T> dual = getDual();
			return "(" + field.power(dual.scale, 2) + "*X + " + dual.xOffset + ", " + field.power(dual.scale, 3)
					+ "*Y + " + field.multiply(dual.scale, dual.scale, dual.xyOffset) + "X + " + dual.yOffset + ")";
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

		public T getScale() {
			return scale;
		}

		public T getXOffset() {
			return xOffset;
		}

		public T getXYOffset() {
			return xyOffset;
		}

		public T getYOffset() {
			return yOffset;
		}

		@Override
		public ProjectivePoint<T> evaluate(ProjectivePoint<T> point) {
			if (point.getCoord(3).equals(field.zero())) {
				return point;
			}
			T x = point.getDehomogenisedCoord(1, 3);
			T y = point.getDehomogenisedCoord(2, 3);
			T xPrime = field.divide(field.subtract(x, xOffset), field.power(scale, 2));
			T yPrime = field.divide(
					field.subtract(y, field.add(field.multiply(xyOffset, field.power(scale, 2), xPrime), yOffset)),
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
			T inverseXYOffset = field.multiply(-1, xyOffset, field.power(inverseScale, 3));
			T inverseYOffset = field.multiply(field.add(field.multiply(xOffset, xyOffset), field.multiply(-1, yOffset)),
					field.power(inverseScale, 3));
			return new Isomorphism<>(range, inverseScale, inverseXOffset, inverseXYOffset, inverseYOffset);
		}

		@Override
		public List<ProjectivePoint<T>> kernelGenerators() {
			return Collections.emptyList();
		}
	}

	public static <T extends Element<T>> Isomorphism<T> concatIsomorphisms(Isomorphism<T> first,
			Isomorphism<T> second) {
		Field<T> field = first.field;
		if (!first.getRange().equals(second.getDomain())) {
			throw new ArithmeticException("Morphisms don't concatenate");
		}
		Isomorphism<T> result = new Isomorphism<>(first.domain, field.multiply(first.scale, second.scale),
				field.add(field.multiply(first.scale, first.scale, second.xOffset), first.xOffset),
				field.add(field.multiply(first.scale, second.xyOffset), first.xyOffset),
				field.add(field.multiply(field.power(first.scale, 3), second.yOffset),
						field.multiply(field.power(first.scale, 2), first.xyOffset, second.xOffset), first.yOffset));
		if (!result.getRange().equals(second.getRange())) {
			throw new ArithmeticException("Concatenation formula wrong!");
		}
		return result;
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
	private Map<Integer, Polynomial<T>> divisionPolynomials;
	private Polynomial<T> definingPolynomial;
	private BigInteger numberOfPoints;
	private boolean superSingular;
	private boolean superSingularDeterminied;
	private Map<IntE, List<ProjectivePoint<T>>> torsionPointBasis;
	private Isomorphism<T> weierstrassForm;
	private GenericProjectiveScheme<T> asGenericProjectiveScheme;
	private ProjectiveMorphism<T> xCover;
	private int minimalDefinedDegree = -1;

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
			return new EllipticCurve<>(field, field.zero(), field.one(), field.zero(), field.zero(),
					field.negative(field.inverse(j)));
		}
		Pair<T, T> coefficients = fromJInvariantCoefficients(field, j);
		return new EllipticCurve<>(field, coefficients.getFirst(), coefficients.getSecond());
	}

	public static <T extends Element<T>, S extends AlgebraicExtensionElement<T, S>, FE extends FieldExtension<T, S, FE>> EllipticCurve<S> extendBaseField(
			FE extension, EllipticCurve<T> curve) {
		return new EllipticCurve<>(extension, extension.getEmbedding(curve.getA1()),
				extension.getEmbedding(curve.getA2()), extension.getEmbedding(curve.getA3()),
				extension.getEmbedding(curve.getA4()), extension.getEmbedding(curve.getA6()));
	}

	public static <T extends Element<T>, S extends AlgebraicExtensionElement<T, S>, FE extends FieldExtension<T, S, FE>> EllipticCurve<S> extendBaseField(
			FieldEmbedding<T, S, FE> extension, EllipticCurve<S> curve) {
		return new EllipticCurve<>(extension.getField(), extension.getEmbedding(curve.getA1()),
				extension.getEmbedding(curve.getA2()), extension.getEmbedding(curve.getA3()),
				extension.getEmbedding(curve.getA4()), extension.getEmbedding(curve.getA6()));
	}

	public static <T extends Element<T>> EllipticCurve<T> fromLegendre(Field<T> field, T lambda) {
		return new EllipticCurve<>(field, field.zero(), field.negative(field.add(field.one(), lambda)), field.zero(),
				lambda, field.zero());
	}

	public static <T extends Element<T>> EllipticCurve<T> fromJInvariantInLegendgeForm(Field<T> field, T j) {
		UnivariatePolynomialRing<T> polynomials = field.getUnivariatePolynomialRing();
		UnivariatePolynomial<T> lhs = polynomials.multiply(j,
				polynomials.power(polynomials.subtract(polynomials.getVarPower(2), polynomials.getVar()), 2));
		UnivariatePolynomial<T> rhs = polynomials.toUnivariate(polynomials.multiply(256,
				polynomials.power(polynomials
						.add(polynomials.subtract(polynomials.getVarPower(2), polynomials.getVar()), polynomials.one()),
						3)));
		Polynomial<T> toSolve = polynomials.subtract(rhs, lhs);
		Map<T, Integer> roots = field.roots(toSolve);
		if (roots.isEmpty()) {
			throw new ArithmeticException("Could not find any roots");
		}
		return fromLegendre(field, roots.keySet().iterator().next());
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
		NFE tauNumerator = field.maximalOrder().getNumerator(tau);
		NFE tauDenominator = field.maximalOrder().getDenominator(tau);
		IntE discriminant = q.subtract(
				q.multiply(field.trace(field.multiply(tauNumerator, tauNumerator)),
						field.trace(field.multiply(tauDenominator, tauDenominator))),
				q.power(field.trace(field.multiply(tauNumerator, tauDenominator)), 2)).asInteger();
		IntE conductor = z.divideChecked(discriminant, field.discriminant());
		NFE orderGenerator = field.multiply(conductor, field.maximalOrder().getModuleGenerators().get(1));
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
		EllipticCurve<NFE> curve = fromJInvariant(classField, jInvariant);
		IntE denominator = z.one();
		denominator = z.lcm(denominator(classField.maximalOrder().getIntegerDenominator(curve.getA1()), 1),
				denominator);
		denominator = z.lcm(denominator(classField.maximalOrder().getIntegerDenominator(curve.getA2()), 2),
				denominator);
		denominator = z.lcm(denominator(classField.maximalOrder().getIntegerDenominator(curve.getA3()), 3),
				denominator);
		denominator = z.lcm(denominator(classField.maximalOrder().getIntegerDenominator(curve.getA4()), 4),
				denominator);
		denominator = z.lcm(denominator(classField.maximalOrder().getIntegerDenominator(curve.getA6()), 6),
				denominator);
		curve = new Isomorphism<>(curve, classField.inverse(classField.getEmbedding(denominator)), classField.zero(),
				classField.zero(), classField.zero()).getRange();
		curve = curve.minimalModel(classField).get().getRange();
		// C/M -> C/M, phi(z) = a*z
		// phi(M) subset M
		// <a, a*tau> subset M
		EllipticCurve<NFE> extendedCurve = rational ? extendBaseField(rayClassField, curve) : curve;
		Isogeny<NFE> endomorphism = new Isogeny<NumberField.NFE>() {

			@Override
			public int compareTo(Isogeny<NFE> o) {
				// TODO Auto-generated method stub
				return 0;
			}

			@Override
			public List<ProjectivePoint<NFE>> kernelGenerators() {
				// TODO Auto-generated method stub
				return null;
			}

			@Override
			public EllipticCurve<NFE> getRange() {
				return extendedCurve;
			}

			@Override
			public Isogeny<NFE> getDual() {
				// TODO Auto-generated method stub
				return null;
			}

			@Override
			public EllipticCurve<NFE> getDomain() {
				return extendedCurve;
			}

			@Override
			public BigInteger getDegree() {
				// TODO Auto-generated method stub
				return null;
			}

			@Override
			public ProjectivePoint<NFE> evaluate(ProjectivePoint<NFE> point) {
				// TODO Auto-generated method stub
				return null;
			}

			@Override
			public ProjectiveMorphism<NFE> asMorphism() {
				// TODO Auto-generated method stub
				return null;
			}
		};
		return new WithComplexMultiplication(curve, null, tau, field, rayClassField);
	}

	private static IntE denominator(IntE a, int weight) {
		Integers z = Integers.z();
		FactorizationResult<IntE, IntE> factors = z.uniqueFactorization(a);
		IntE result = z.one();
		for (IntE prime : factors.primeFactors()) {
			result = z.multiply(z.power(prime, MiscAlgorithms.DivRoundUp(factors.multiplicity(prime), weight)), result);
		}
		return result;
	}

	public static enum ReductionType {
		GOOD_REDUCTION, POTENTIAL_GOOD_REDUCTION, BAD_REDUCTION;
	}

	public static enum KodairaSymbol {
		I0(ReductionType.GOOD_REDUCTION), I1(ReductionType.BAD_REDUCTION), I2(ReductionType.BAD_REDUCTION),
		Iv(ReductionType.BAD_REDUCTION), mIv(ReductionType.BAD_REDUCTION), II(ReductionType.POTENTIAL_GOOD_REDUCTION),
		III(ReductionType.POTENTIAL_GOOD_REDUCTION), IV(ReductionType.POTENTIAL_GOOD_REDUCTION),
		I0star(ReductionType.POTENTIAL_GOOD_REDUCTION), Ivstar(ReductionType.POTENTIAL_GOOD_REDUCTION),
		IVstar(ReductionType.POTENTIAL_GOOD_REDUCTION), IIIstar(ReductionType.POTENTIAL_GOOD_REDUCTION),
		IIstar(ReductionType.POTENTIAL_GOOD_REDUCTION);

		private KodairaSymbol(ReductionType type) {
			this.reductionType = type;
		}

		private final ReductionType reductionType;

		public ReductionType reductionType() {
			return reductionType;
		}
	}

	public static <T extends Element<T>> EllipticCurve<T> fromPoints(Field<T> field, List<ProjectivePoint<T>> points) {
		return fromPoints(field, points, Optional.empty(), Optional.empty(), Optional.empty(), Optional.empty(),
				Optional.empty());
	}

	public static <T extends Element<T>> EllipticCurve<T> fromWeierstrassPoints(Field<T> field,
			List<ProjectivePoint<T>> points) {
		if (field.characteristic().compareTo(BigInteger.valueOf(3)) <= 0) {
			throw new ArithmeticException("Not possible for small characteristic!");
		}
		return fromPoints(field, points, Optional.of(field.zero()), Optional.of(field.zero()),
				Optional.of(field.zero()), Optional.empty(), Optional.empty());
	}

	public static <T extends Element<T>> EllipticCurve<T> fromPoints(Field<T> field, List<ProjectivePoint<T>> points,
			Optional<T> a1, Optional<T> a2, Optional<T> a3, Optional<T> a4, Optional<T> a6) {
		List<List<T>> equations = new ArrayList<>();
		List<T> rhs = new ArrayList<>();
		for (ProjectivePoint<T> point : points) {
			if (point.equals(new ProjectivePoint<>(field, field.zero(), field.one(), field.zero()))) {
				continue;
			}
			List<T> row = new ArrayList<>();
			T x = point.getDehomogenisedCoord(1, 3);
			T y = point.getDehomogenisedCoord(2, 3);
			T target = field.subtract(field.power(y, 2), field.power(x, 3));
			if (a1.isPresent()) {
				target = field.add(target, field.multiply(a1.get(), x, y));
			} else {
				row.add(field.negative(field.multiply(x, y)));
			}
			if (a2.isPresent()) {
				target = field.subtract(target, field.multiply(a2.get(), x, x));
			} else {
				row.add(field.multiply(x, x));
			}
			if (a3.isPresent()) {
				target = field.add(target, field.multiply(a3.get(), y));
			} else {
				row.add(field.negative(y));
			}
			if (a4.isPresent()) {
				target = field.subtract(target, field.multiply(a4.get(), x));
			} else {
				row.add(x);
			}
			if (a6.isPresent()) {
				target = field.subtract(target, a6.get());
			} else {
				row.add(field.one());
			}
			equations.add(row);
			rhs.add(target);
		}
		Matrix<T> asMatrix = new Matrix<>(equations);
		Vector<T> asVector = new Vector<>(rhs);
		if (!field.isSubModuleMember(asMatrix, asVector)) {
			throw new ArithmeticException("Points do not share an elliptic curve with the given constraints!");
		}
		if (!field.syzygyProblem(asMatrix).isEmpty()) {
			throw new ArithmeticException("Points do not uniquely identify an elliptic curve!");
		}
		Vector<T> solution = field.asSubModuleMember(asMatrix, asVector);
		int index = 1;
		T actualA1;
		if (a1.isEmpty()) {
			actualA1 = solution.get(index);
			index++;
		} else {
			actualA1 = a1.get();
		}
		T actualA2;
		if (a2.isEmpty()) {
			actualA2 = solution.get(index);
			index++;
		} else {
			actualA2 = a2.get();
		}
		T actualA3;
		if (a3.isEmpty()) {
			actualA3 = solution.get(index);
			index++;
		} else {
			actualA3 = a3.get();
		}
		T actualA4;
		if (a4.isEmpty()) {
			actualA4 = solution.get(index);
			index++;
		} else {
			actualA4 = a4.get();
		}
		T actualA6;
		if (a6.isEmpty()) {
			actualA6 = solution.get(index);
			index++;
		} else {
			actualA6 = a6.get();
		}
		return new EllipticCurve<>(field, actualA1, actualA2, actualA3, actualA4, actualA6);
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
		this.projectiveRing = AbstractPolynomialRing.getPolynomialRing(field, 3, Monomial.GREVLEX);
		this.affineRing = AbstractPolynomialRing.getPolynomialRing(this.field, 2, Monomial.REVLEX);
		PolynomialRing<T> r = this.affineRing;
		Polynomial<T> a6p = r.getEmbedding(a6);
		Polynomial<T> a4x = r.getEmbedding(a4, new int[] { 1, 0 });
		Polynomial<T> a2x2 = r.getEmbedding(a2, new int[] { 2, 0 });
		Polynomial<T> x3 = r.getVarPower(1, 3);
		Polynomial<T> y2 = r.getVarPower(2, 2);
		Polynomial<T> a1xy = r.getEmbedding(a1, new int[] { 1, 1 });
		Polynomial<T> a3y = r.getEmbedding(a3, new int[] { 0, 1 });
		Polynomial<T> rhs = r.add(x3, a2x2, a4x, a6p);
		Polynomial<T> lhs = r.add(y2, a1xy, a3y);
		this.definingPolynomial = r.subtract(rhs, lhs);
		this.univariateRing = this.field.getUnivariatePolynomialRing();
		this.torsionPointBasis = new TreeMap<>();
	}

	public GenericProjectiveScheme<T> asGenericProjectiveScheme() {
		if (asGenericProjectiveScheme == null) {
			PolynomialRing<T> ring = projectiveRing;
			Polynomial<T> f = ring.getEmbedding(field.one(), new int[] { 3, 0, 0 });
			f = ring.add(f, ring.getEmbedding(field.negative(field.one()), new int[] { 0, 2, 1 }));
			f = ring.add(f, ring.getEmbedding(field.negative(a1), new int[] { 1, 1, 1 }));
			f = ring.add(f, ring.getEmbedding(field.negative(a3), new int[] { 0, 1, 2 }));
			f = ring.add(f, ring.getEmbedding(a2, new int[] { 2, 0, 1 }));
			f = ring.add(f, ring.getEmbedding(a4, new int[] { 1, 0, 2 }));
			f = ring.add(f, ring.getEmbedding(a6, new int[] { 0, 0, 3 }));
			asGenericProjectiveScheme = new GenericProjectiveScheme<>(field, ring, Collections.singletonList(f));
		}
		return asGenericProjectiveScheme;
	}

	public T jInvariant() {
		return j;
	}

	public T discriminant() {
		return discriminant;
	}

	public Isomorphism<T> identity() {
		return new Isomorphism<>(this, field.one(), field.zero(), field.zero(), field.zero());
	}

	@Override
	public ProjectiveMorphism<T> identityMorphism() {
		return identity().asMorphism();
	}

	private List<Polynomial<T>> linearTransforms(PolynomialRing<T> r) {
		List<Polynomial<T>> equations = new ArrayList<>();
		Polynomial<T> s = r.getVar(1);
		Polynomial<T> x = r.getVar(2);
		Polynomial<T> xy = r.getVar(3);
		Polynomial<T> y = r.getVar(4);
		Polynomial<T> a1 = r.getEmbedding(getA1());
		Polynomial<T> a2 = r.getEmbedding(getA2());
		Polynomial<T> a3 = r.getEmbedding(getA3());
		Polynomial<T> a4 = r.getEmbedding(getA4());
		Polynomial<T> a6 = r.getEmbedding(getA6());
		equations.add(r.multiply(s, r.add(r.multiply(2, xy), a1)));
		equations.add(
				r.multiply(r.power(s, 2), r.add(a2, r.multiply(-1, xy, a1), r.multiply(3, x), r.multiply(-1, xy, xy))));
		equations.add(r.multiply(r.power(s, 3), r.add(a3, r.multiply(x, a1), r.multiply(2, y))));
		equations.add(r.multiply(r.power(s, 4), r.add(a4, r.multiply(-1, xy, a3), r.multiply(2, x, a2),
				r.add(r.multiply(-1, r.add(y, r.multiply(x, xy)), a1), r.multiply(3, x, x), r.multiply(-2, xy, y)))));
		equations.add(r.multiply(r.power(s, 6), r.add(r.add(a6, r.multiply(x, a4), r.multiply(x, x, a2)),
				r.add(r.power(x, 3), r.multiply(-1, y, a3), r.multiply(-1, y, y), r.multiply(-1, x, y, a1)))));
		return equations;
	}

	public List<Isomorphism<T>> getIsomorphisms(EllipticCurve<T> range) {
		if (!this.jInvariant().equals(range.jInvariant())) {
			return Collections.emptyList();
		}
//		a1' = 2*s^3*xyO + a1*s
//		a2' = -1*s^6*xyO^2 + -1*a1*s^4*xyO + a2*s^2 + 3*s^2*xO
//		a3' = a1*s^3*xO + a3*s^3 + 2*s^3yO
//		a4' = -1*a1s^6*xO*xyO + -1*a3*s^6*xyO + -2*s^6*xyO*yO + 2*a2*s^4*xO + 3*s^4*xO^2 + -1*a1*s^4*yO + a4*s^4
//		a6' = a2*s^6*xO^2 + s^6*xO^3 + -1*a1*s^6*xO*yO + a4*s^6*xO + -1*a3*s^6*yO + -1*s^6*yO^2 + a6*s^6
		PolynomialRing<T> r = AbstractPolynomialRing.getPolynomialRing(field, Monomial.REVLEX,
				new String[] { "u", "r", "s", "t" });
		List<Polynomial<T>> equations = linearTransforms(r);
		equations.set(0, r.subtract(equations.get(0), r.getEmbedding(a1)));
		equations.set(1, r.subtract(equations.get(1), r.getEmbedding(a2)));
		equations.set(2, r.subtract(equations.get(2), r.getEmbedding(a3)));
		equations.set(3, r.subtract(equations.get(3), r.getEmbedding(a4)));
		equations.set(4, r.subtract(equations.get(4), r.getEmbedding(a6)));
		List<AffinePoint<T>> solved = r.solve(equations);
		List<Isomorphism<T>> isomorphisms = new ArrayList<>();
		for (AffinePoint<T> root : solved) {
			isomorphisms.add(new Isomorphism<>(this, field.inverse(root.getCoord(1)), root.getCoord(2),
					root.getCoord(3), root.getCoord(4)));
		}
		return isomorphisms;
	}

	public List<Isomorphism<T>> getAutomorphisms() {
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
		List<Isomorphism<T>> isomorphisms = frobenius.getRange().getIsomorphisms(this);
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

	public T getB2() {
		return b2;
	}

	public T getB4() {
		return b4;
	}

	public T getB6() {
		return b6;
	}

	public T getB8() {
		return b8;
	}

	public T getC4() {
		return c4;
	}

	public T getC6() {
		return c6;
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

	public Isomorphism<T> getWeierstrassForm() {
		if (weierstrassForm != null) {
			return weierstrassForm;
		}
		T zero = field.zero();
		T one = field.one();
		if (field.characteristic().equals(BigInteger.TWO)) {
			if (!a1.equals(zero)) {
				if (a1.equals(one) && a3.equals(zero) && a4.equals(zero)) {
					weierstrassForm = this.identity();
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
			if (a1.equals(zero) && a3.equals(zero) && (a2.equals(zero) || a4.equals(zero))) {
				weierstrassForm = this.identity();
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
			if (!transform.isEmpty()) {
				weierstrassForm = new Isomorphism<>(this, one, transform.get(0).getCoord(1),
						transform.get(0).getCoord(2), transform.get(0).getCoord(3));
			} else {
				requiredEquations.clear();
				requiredEquations.add(r3.getEmbedding(r4.partiallyEvaluate(equations.get(0), eval), map));
				requiredEquations.add(r3.getEmbedding(r4.partiallyEvaluate(equations.get(2), eval), map));
				requiredEquations.add(r3.getEmbedding(r4.partiallyEvaluate(equations.get(3), eval), map));
				transform = r3.solve(requiredEquations);
				weierstrassForm = new Isomorphism<>(this, one, transform.get(0).getCoord(1),
						transform.get(0).getCoord(2), transform.get(0).getCoord(3));

			}
			return weierstrassForm;
		}
		if (a1.equals(zero) && a3.equals(zero) && a2.equals(zero)) {
			weierstrassForm = this.identity();
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

	private static <T extends Element<T>> List<T> otherLambdas(Field<T> field, T lambda) {
		List<T> result = new ArrayList<>();
		result.add(lambda);
		result.add(field.subtract(field.one(), lambda));
		result.add(field.inverse(lambda));
		result.add(field.inverse(field.subtract(field.one(), lambda)));
		result.add(field.subtract(field.one(), field.inverse(lambda)));
		result.add(field.inverse(field.subtract(field.one(), field.inverse(lambda))));
		return result;
	}

	public static class LegendreFormIsomorphism<T extends Element<T>> {
		private Isomorphism<T> isomorphism;
		private T lambda;
		private List<T> otherLambdas;

		private LegendreFormIsomorphism(Isomorphism<T> isomorphism, T lambda) {
			this.isomorphism = isomorphism;
			this.lambda = lambda;
			this.otherLambdas = otherLambdas(isomorphism.field, lambda);
		}

		public Isomorphism<T> getIsomorphism() {
			return isomorphism;
		}

		public T getLambda() {
			return lambda;
		}

		public List<T> getOtherLambdas() {
			return otherLambdas;
		}
	}

	public Optional<LegendreFormIsomorphism<T>> getLengendreFormIsomorphism() {
		if (field.characteristic().equals(BigInteger.TWO)) {
			throw new ArithmeticException("Characteristic equals 2!");
		}
		if (!getA1().equals(field.zero()) || !getA3().equals(field.zero())) {
			Isomorphism<T> weierstrass = getWeierstrassForm();
			Optional<LegendreFormIsomorphism<T>> legendre = weierstrass.getRange().getLengendreFormIsomorphism();
			return legendre.map((LegendreFormIsomorphism<T> form) -> {
				return new LegendreFormIsomorphism<>(concatIsomorphisms(weierstrass, form.getIsomorphism()),
						form.getLambda());
			});
		}
		UnivariatePolynomial<T> polynomial = getAdjustedRhs();
		Map<T, Integer> roots = field.roots(polynomial);
		if (roots.size() != 3) {
			return Optional.empty();
		}
		Iterator<T> it = roots.keySet().iterator();
		T root1 = it.next();
		T root2 = it.next();
		T root3 = it.next();
		T scaleSqr = field.subtract(root2, root1);
		if (!field.hasSqrt(scaleSqr)) {
			return Optional.empty();
		}
		T scale = field.sqrt(scaleSqr).keySet().iterator().next();
		T xOffset = root1;
		Isomorphism<T> isomorphism = new Isomorphism<>(this, scale, xOffset, field.zero(), field.zero());
		T lambda = getLambda(field, root1, root2, root3);
		if (!isomorphism.getRange().jInvariant().equals(fromLegendre(field, lambda).jInvariant())) {
			throw new ArithmeticException("Calculation wrong!");
		}
		return Optional.of(new LegendreFormIsomorphism<>(isomorphism, lambda));
	}

	public static class ExtensionLegendreFormIsomorphism<T extends Element<T>, B extends Element<B>, E extends AlgebraicExtensionElement<B, E>, FE extends FieldExtension<B, E, FE>> {
		private Extension<T, B, E, FE> extension;
		private E lambda;
		private List<E> otherLambdas;
		private EllipticCurve<E> curve;
		private MathMap<ProjectivePoint<T>, ProjectivePoint<E>> embedding;

		private ExtensionLegendreFormIsomorphism(Extension<T, B, E, FE> extension, EllipticCurve<E> curve, E lambda,
				MathMap<ProjectivePoint<T>, ProjectivePoint<E>> embedding) {
			this.extension = extension;
			this.lambda = lambda;
			this.curve = curve;
			this.embedding = embedding;
			this.otherLambdas = otherLambdas(extension.extension(), lambda);
		}

		public EllipticCurve<E> getCurve() {
			return curve;
		}

		public E getLambda() {
			return lambda;
		}

		public Extension<T, B, E, FE> getExtension() {
			return extension;
		}

		public MathMap<ProjectivePoint<T>, ProjectivePoint<E>> getEmbedding() {
			return embedding;
		}

		public List<E> getOtherLambdas() {
			return otherLambdas;
		}
	}

	public ExtensionLegendreFormIsomorphism<T, ?, ?, ?> getExtensionLengendreFormIsomorphism() {
		return getExtensionLengendreFormIsomorphism(field.getExtension(field.getUnivariatePolynomialRing().getVar()));
	}

	public <B extends Element<B>, E extends AlgebraicExtensionElement<B, E>, FE extends FieldExtension<B, E, FE>> ExtensionLegendreFormIsomorphism<T, B, E, FE> getExtensionLengendreFormIsomorphism(
			Extension<T, B, E, FE> trivialExtension) {
		if (field.characteristic().equals(BigInteger.TWO)) {
			throw new ArithmeticException("Characteristic equals 2!");
		}
		if (!getA1().equals(field.zero()) || !getA3().equals(field.zero())) {
			Isomorphism<T> weierstrass = getWeierstrassForm();
			ExtensionLegendreFormIsomorphism<T, B, E, FE> legendre = weierstrass.getRange()
					.getExtensionLengendreFormIsomorphism(trivialExtension);
			return new ExtensionLegendreFormIsomorphism<T, B, E, FE>(legendre.getExtension(), legendre.getCurve(),
					legendre.getLambda(), new FunctionMathMap<>((ProjectivePoint<T> point) -> {
						return legendre.getEmbedding().evaluate(weierstrass.evaluate(point));
					}));
		}
		FE extension = trivialExtension.extension();
		UnivariatePolynomial<E> polynomial = extension.getUnivariatePolynomialRing().getEmbedding(getAdjustedRhs(),
				trivialExtension.embeddingMap());
		SplittingFieldResult<B, E, FE> splittingField = extension.getSplittingField(polynomial);
		FieldEmbedding<B, E, FE> embeddedExtension = splittingField.getExtension();
		extension = embeddedExtension.getField();
		E root1 = splittingField.getConjugates().get(0);
		E root2 = splittingField.getConjugates().get(1);
		E root3 = splittingField.getConjugates().get(2);
		E scaleSqr = extension.subtract(root2, root1);
		E scale;
		E xOffset = root1;
		E lambda = getLambda(extension, root1, root2, root3);
		if (!extension.hasSqrt(scaleSqr)) {
			List<E> coeffs = new ArrayList<>();
			coeffs.add(extension.negative(scaleSqr));
			coeffs.add(extension.zero());
			coeffs.add(extension.one());
			UnivariatePolynomial<E> embedded = extension.getUnivariatePolynomialRing().getPolynomial(coeffs);
			FieldEmbedding<B, E, FE> sqrtExtension = embeddedExtension.getField().getEmbeddedExtension(embedded);
			embeddedExtension = new FieldEmbedding<>(embeddedExtension, sqrtExtension);
			extension = embeddedExtension.getField();
			scale = sqrtExtension.getGenerator();
			xOffset = sqrtExtension.getEmbedding(xOffset);
			lambda = sqrtExtension.getEmbedding(lambda);
		} else {
			scale = extension.sqrt(scaleSqr).keySet().iterator().next();
		}
		Extension<T, B, E, FE> lambdaExtension = trivialExtension.extendFurther(embeddedExtension);
		EllipticCurveExtensionResult<T, B, E, FE> ext = extendBaseField(lambdaExtension);
		Isomorphism<E> isomorphism = new Isomorphism<E>(ext.getCurve(), scale, xOffset, extension.zero(),
				extension.zero());
		return new ExtensionLegendreFormIsomorphism<>(lambdaExtension, isomorphism.getRange(), lambda,
				new FunctionMathMap<>((ProjectivePoint<T> point) -> {
					return isomorphism.evaluate(ext.getEmbedding().evaluate(point));
				}));
	}

	public static class LegendreForm<T extends Element<T>> {
		private EllipticCurve<T> curve;
		private T lambda;
		private List<T> otherLambdas;

		private LegendreForm(EllipticCurve<T> curve, T lambda) {
			this.curve = curve;
			this.lambda = lambda;
			this.otherLambdas = otherLambdas(curve.getField(), lambda);
		}

		public EllipticCurve<T> getCurve() {
			return curve;
		}

		public T getLambda() {
			return lambda;
		}

		public List<T> getOtherLambdas() {
			return otherLambdas;
		}
	}

	public Optional<LegendreForm<T>> getLengendreForm() {
		if (field.characteristic().equals(BigInteger.TWO)) {
			throw new ArithmeticException("Characteristic equals 2!");
		}
		if (!getA1().equals(field.zero()) || !getA3().equals(field.zero())) {
			return fromJInvariant(field, jInvariant()).getLengendreForm();
		}
		List<T> rhs = new ArrayList<>();
		rhs.add(getA6());
		rhs.add(getA4());
		rhs.add(getA2());
		rhs.add(field.one());
		UnivariatePolynomialRing<T> polynomialRing = field.getUnivariatePolynomialRing();
		UnivariatePolynomial<T> polynomial = polynomialRing.getPolynomial(rhs);
		Map<T, Integer> roots = field.roots(polynomial);
		if (roots.size() != 3) {
			return Optional.empty();
		}
		Iterator<T> it = roots.keySet().iterator();
		T root1 = it.next();
		T root2 = it.next();
		T root3 = it.next();
		T lambda = getLambda(field, root1, root2, root3);
		EllipticCurve<T> curve = fromLegendre(field, lambda);
		if (!jInvariant().equals(curve.jInvariant())) {
			throw new ArithmeticException("Calculation was wrong!");
		}
		return Optional.of(new LegendreForm<>(curve, lambda));
	}

	public static class ExtensionLegendreForm<T extends Element<T>, B extends Element<B>, E extends AlgebraicExtensionElement<B, E>, FE extends FieldExtension<B, E, FE>> {
		private Extension<T, B, E, FE> extension;
		private E lambda;
		private List<E> otherLambdas;
		private EllipticCurve<E> curve;

		private ExtensionLegendreForm(Extension<T, B, E, FE> extension, EllipticCurve<E> curve, E lambda) {
			this.extension = extension;
			this.lambda = lambda;
			this.otherLambdas = otherLambdas(extension.extension(), lambda);
			this.curve = curve;
		}

		public EllipticCurve<E> getCurve() {
			return curve;
		}

		public E getLambda() {
			return lambda;
		}

		public Extension<T, B, E, FE> getExtension() {
			return extension;
		}

		public List<E> getOtherLambdas() {
			return otherLambdas;
		}
	}

	public ExtensionLegendreForm<T, ?, ?, ?> getExtensionLengendreForm() {
		return getExtensionLengendreForm(field.getExtension(field.getUnivariatePolynomialRing().getVar()));
	}

	public <B extends Element<B>, E extends AlgebraicExtensionElement<B, E>, FE extends FieldExtension<B, E, FE>> ExtensionLegendreForm<T, B, E, FE> getExtensionLengendreForm(
			Extension<T, B, E, FE> trivialExtension) {
		if (field.characteristic().equals(BigInteger.TWO)) {
			throw new ArithmeticException("Characteristic equals 2!");
		}
		if (!getA1().equals(field.zero()) || !getA3().equals(field.zero())) {
			return fromJInvariant(field, jInvariant()).getExtensionLengendreForm(trivialExtension);
		}
		FE extension = trivialExtension.extension();
		UnivariatePolynomial<E> polynomial = extension.getUnivariatePolynomialRing().getEmbedding(getAdjustedRhs(),
				trivialExtension.embeddingMap());
		SplittingFieldResult<B, E, FE> splittingField = extension.getSplittingField(polynomial);
		FieldEmbedding<B, E, FE> embeddedExtension = splittingField.getExtension();
		Extension<T, B, E, FE> lambdaExtension = trivialExtension.extendFurther(embeddedExtension);
		E lambda = getLambda(embeddedExtension.getField(), splittingField.getConjugates().get(0),
				splittingField.getConjugates().get(1), splittingField.getConjugates().get(2));
		return new ExtensionLegendreForm<>(lambdaExtension, fromLegendre(embeddedExtension.getField(), lambda), lambda);
	}

	private static <T extends Element<T>> T getLambda(Field<T> field, T root1, T root2, T root3) {
		T scaleSqr = field.subtract(root2, root1);
		T xOffset = root1;
		return field.divide(field.subtract(root3, xOffset), scaleSqr);
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
			EllipticCurve<T> twistedCurve = new EllipticCurve<>(field, a1, field.add(a2, field.multiply(twist, a1, a1)),
					a3, a4, field.add(a6, field.multiply(twist, a3, a3)));
			if (numberOfPoints != null) {
				twistedCurve.numberOfPoints = field.getNumberOfElements().add(BigInteger.ONE).multiply(BigInteger.TWO)
						.subtract(numberOfPoints);
			}
			return twistedCurve;
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
			Isogeny<T> isogeny = new Isomorphism<>(this, field.one(), field.zero(), transform.get(0).getCoord(1),
					transform.get(0).getCoord(2));
			setNumberOfPointsFrom(this, isogeny.getDual());
			return isogeny.getRange().getQuadraticTwist();

		}
		T twist;
		do {
			twist = field.getRandomElement();
		} while (field.hasSqrt(twist));
		EllipticCurve<T> twistedCurve = new EllipticCurve<>(field, field.zero(), field.multiply(twist, a2),
				field.zero(), field.multiply(twist, twist, a4), field.multiply(field.power(twist, 3), a6));
		if (numberOfPoints != null) {
			twistedCurve.numberOfPoints = field.getNumberOfElements().add(BigInteger.ONE).multiply(BigInteger.TWO)
					.subtract(numberOfPoints);
		}
		return twistedCurve;
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

//	private T finiteFieldNorm(T t) {
//		return finiteFieldNorm(t, 1);
//	}
//
//	private T finiteFieldNorm(T t, int degree) {
//		BigInteger q = field.getNumberOfElements();
//		BigInteger p = field.characteristic().pow(degree);
//		T result = field.one();
//		T power = t;
//		while (!q.equals(BigInteger.ONE)) {
//			result = field.multiply(power, result);
//			power = field.power(power, p);
//			q = q.divide(p);
//		}
//		return result;
//	}

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
			if (field.add(field.multiply(2, py), field.multiply(a1, px), a3).equals(field.zero()) || !py.equals(qy)) {
				return field.subtract(atx, px);
			}
			s = field.divide(
					field.add(field.multiply(3, px, px), field.multiply(2, a2, px), a4, field.multiply(-1, a1, py)),
					field.add(field.multiply(2, py), field.multiply(a1, px), a3));
		} else {
			s = field.divide(field.subtract(py, qy), field.subtract(px, qx));
		}
		T num = field.subtract(aty, field.add(py, field.multiply(s, field.subtract(atx, px))));
		T denom = field.subtract(field.add(atx, px, qx, a2), field.add(field.multiply(s, s), field.multiply(s, a1)));
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

	private T weilPairingExtension(BigInteger n, ProjectivePoint<T> p, ProjectivePoint<T> q) {
		UnivariatePolynomial<T> minimalPolynomial;
		if (field.characteristic().equals(BigInteger.TWO)) {
			T nonZeroTrace;
			do {
				nonZeroTrace = field.getRandomElement();
			} while (finiteFieldTrace(nonZeroTrace).equals(field.zero()));
			List<T> coeffs = new ArrayList<>();
			coeffs.add(nonZeroTrace);
			coeffs.add(field.one());
			coeffs.add(field.one());
			minimalPolynomial = univariateRing.getPolynomial(coeffs);
		} else {
			T nonSquare;
			do {
				nonSquare = field.getRandomElement();
			} while (field.hasSqrt(nonSquare));
			List<T> coeffs = new ArrayList<>();
			coeffs.add(field.negative(nonSquare));
			coeffs.add(field.zero());
			coeffs.add(field.one());
			minimalPolynomial = univariateRing.getPolynomial(coeffs);
		}
		return weilPairingExtension(n, p, q, field.getExtension(minimalPolynomial));
	}

	private <B extends Element<B>, E extends AlgebraicExtensionElement<B, E>, FE extends FieldExtension<B, E, FE>> T weilPairingExtension(
			BigInteger n, ProjectivePoint<T> p, ProjectivePoint<T> q, Extension<T, B, E, FE> extension) {
		EllipticCurveExtensionResult<T, B, E, FE> extendedCurve = extendBaseField(extension);
		ProjectivePoint<E> extendedP = extendedCurve.getEmbedding().evaluate(p);
		ProjectivePoint<E> extendedQ = extendedCurve.getEmbedding().evaluate(q);
		E pairing = extendedCurve.getCurve().weilPairing(n, extendedP, extendedQ);
		return extension.retractionMap().evaluate(pairing);
	}

	// Silverman XI 8.
	public T weilPairing(BigInteger n, ProjectivePoint<T> p, ProjectivePoint<T> q) {
		if (!multiply(n, p).equals(neutral()) || !multiply(n, q).equals(neutral())) {
			throw new ArithmeticException("Points have wrong order!");
		}
		ProjectivePoint<T> s;
		int tries = 0;
		do {
			s = getRandomElement();
			tries++;
			if (tries > 8) {
				return weilPairingExtension(n, p, q);
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

	public <I extends Element<I>, S extends Element<S>> Isomorphism<T> integerModel(GlobalField<T, I, S> field) {
		DedekindRing<I, T, S> integers = field.ringOfIntegers();
		I denominator = integers.getDenominator(a1);
		denominator = integers.lcm(integers.getDenominator(a2), denominator);
		denominator = integers.lcm(integers.getDenominator(a3), denominator);
		denominator = integers.lcm(integers.getDenominator(a4), denominator);
		denominator = integers.lcm(integers.getDenominator(a6), denominator);
		return new Isomorphism<>(this, field.inverse(field.getInteger(denominator)), field.zero(), field.zero(),
				field.zero());
	}

	public <I extends Element<I>, S extends Element<S>> Optional<Isomorphism<T>> minimalIntegerModel(
			GlobalField<T, I, S> field) {
		Isomorphism<T> integer = integerModel(field);
		return integer.getRange().minimalModel(field)
				.map((Isomorphism<T> minimal) -> concatIsomorphisms(integer, minimal));
	}

	public <I extends Element<I>, S extends Element<S>> Optional<Isomorphism<T>> minimalModel(
			GlobalField<T, I, S> field) {
		if (!field.isInteger(a1) || !field.isInteger(a2) || !field.isInteger(a3) || !field.isInteger(a4)
				|| !field.isInteger(a6)) {
			throw new ArithmeticException("Not an integer curve!");
		}
		DedekindRing<I, T, S> integers = field.ringOfIntegers();
		FactorizationResult<Ideal<I>, Ideal<I>> discriminantFactors = integers
				.idealFactorization(integers.getIdeal(Collections.singletonList(integers.asInteger(discriminant()))));
		Isomorphism<T> minimalModel = identity();
		EllipticCurve<T> curve = this;
		for (Ideal<I> ideal : discriminantFactors.primeFactors()) {
			if (!ideal.isPrincipal()) {
				return Optional.empty();
			}
			DiscreteValuationField<T, S> localField = integers.localize(ideal).localField();
			Isomorphism<T> localMinimalModel = curve.localMinimalModel(localField);
			minimalModel = concatIsomorphisms(minimalModel, localMinimalModel);
			curve = minimalModel.getRange();
			if (!field.isInteger(curve.a1) || !field.isInteger(curve.a2) || !field.isInteger(curve.a3)
					|| !field.isInteger(curve.a4) || !field.isInteger(curve.a6)) {
				System.out.println("broken");
			}
			minimalModel = concatIsomorphisms(minimalModel, curve.simplifyMinimalModel(field));
			curve = minimalModel.getRange();
			if (!field.isInteger(curve.a1) || !field.isInteger(curve.a2) || !field.isInteger(curve.a3)
					|| !field.isInteger(curve.a4) || !field.isInteger(curve.a6)) {
				return Optional.empty();
			}
		}
		return Optional.of(minimalModel);
	}

	private <I extends Element<I>, S extends Element<S>> Isomorphism<T> simplifyMinimalModel(
			GlobalField<T, I, S> field) {
		DedekindRing<I, T, S> integers = field.ringOfIntegers();
		I two = integers.getInteger(2);
		I three = integers.getInteger(3);
		Ideal<I> even = integers.getIdeal(Collections.singletonList(two));
		Ideal<I> multiplesOfThree = integers.getIdeal(Collections.singletonList(three));
		I a1Mod2 = even.residue(integers.asInteger(a1));
		T xyOffset = field.getInteger(integers.divideChecked(integers.subtract(a1Mod2, integers.asInteger(a1)), two));
		T a2Prime = field.add(a2, field.multiply(-1, xyOffset, a1), field.multiply(-1, xyOffset, xyOffset));
		I a2PrimeMod3 = multiplesOfThree.residue(integers.asInteger(a2Prime));
		T xOffset = field
				.getInteger(integers.divideChecked(integers.subtract(a2PrimeMod3, integers.asInteger(a2Prime)), three));
		T a3Prime = field.add(a3, field.multiply(xOffset, a1));
		I a3PrimeMod2 = even.residue(integers.asInteger(a3Prime));
		T yOffset = field
				.getInteger(integers.divideChecked(integers.subtract(a3PrimeMod2, integers.asInteger(a3Prime)), two));
		return new Isomorphism<>(this, field.one(), xOffset, xyOffset, yOffset);
	}

	public <S extends Element<S>> Isomorphism<T> localMinimalModel(DiscreteValuationField<T, S> field) {
		return tatesAlgorithm(field).getTranslation();
	}

	public static class TatesAlgorithmResult<T extends Element<T>, S extends Element<S>> {
		private DiscreteValuationField<T, S> field;
		private Isomorphism<T> translation;
		private KodairaSymbol kodairaSymbol;
		private int value;
		private int multiplicity;
		private ReductionType reductionType;
		private int conductorExponent;
		private int localIndex;
		private EllipticCurve<S> reducedCurve;
		private GenericProjectiveScheme<S> reducedScheme;
		private MathMap<ProjectivePoint<T>, ProjectivePoint<S>> reductionMap;

		private TatesAlgorithmResult(DiscreteValuationField<T, S> field, Isomorphism<T> translation,
				KodairaSymbol kodairaSymbol, int conductorExponent, int localIndex) {
			if (kodairaSymbol == KodairaSymbol.Iv || kodairaSymbol == KodairaSymbol.Ivstar) {
				throw new ArithmeticException("Expected value symbol constructor!");
			}
			this.field = field;
			this.translation = translation;
			this.kodairaSymbol = kodairaSymbol;
			this.value = -1;
			if (kodairaSymbol == KodairaSymbol.I0) {
				value = 0;
			}
			if (kodairaSymbol == KodairaSymbol.I1) {
				value = 1;
			}
			if (kodairaSymbol == KodairaSymbol.I2) {
				value = 2;
			}
			if (kodairaSymbol == KodairaSymbol.I0star) {
				value = 0;
			}
			this.multiplicity = 1;
			this.reductionType = kodairaSymbol.reductionType();
			this.conductorExponent = conductorExponent;
			this.localIndex = localIndex;
			init();
		}

		private TatesAlgorithmResult(DiscreteValuationField<T, S> field, Isomorphism<T> translation,
				KodairaSymbol kodairaSymbol, int value, int conductorExponent, int localIndex) {
			if (kodairaSymbol != KodairaSymbol.Iv && kodairaSymbol != KodairaSymbol.Ivstar) {
				throw new ArithmeticException("Expected non value symbol constructor!");
			}
			if (kodairaSymbol == KodairaSymbol.Iv && value == 0) {
				kodairaSymbol = KodairaSymbol.I0;
			}
			if (kodairaSymbol == KodairaSymbol.Iv && value == 1) {
				kodairaSymbol = KodairaSymbol.I1;
			}
			if (kodairaSymbol == KodairaSymbol.Iv && value == 2) {
				kodairaSymbol = KodairaSymbol.I2;
			}
			if (kodairaSymbol == KodairaSymbol.Ivstar && value == 0) {
				kodairaSymbol = KodairaSymbol.I0star;
			}
			this.field = field;
			this.translation = translation;
			this.kodairaSymbol = kodairaSymbol;
			this.value = value;
			this.multiplicity = 1;
			this.reductionType = kodairaSymbol.reductionType();
			this.conductorExponent = conductorExponent;
			this.localIndex = localIndex;
			init();
		}

		private TatesAlgorithmResult(DiscreteValuationField<T, S> field, Isomorphism<T> translation,
				KodairaSymbol kodairaSymbol, int value, int multiplicity, int conductorExponent, int localIndex) {
			if (kodairaSymbol != KodairaSymbol.mIv) {
				throw new ArithmeticException("Expected multiplicity 1 constructor!");
			}
			if (multiplicity == 1) {
				if (value == 0) {
					kodairaSymbol = KodairaSymbol.I0;
				} else if (value == 1) {
					kodairaSymbol = KodairaSymbol.I1;
				} else if (value == 2) {
					kodairaSymbol = KodairaSymbol.I2;
				} else {
					kodairaSymbol = KodairaSymbol.Iv;
				}
			}
			this.field = field;
			this.translation = translation;
			this.kodairaSymbol = kodairaSymbol;
			this.value = value;
			this.multiplicity = multiplicity;
			this.reductionType = kodairaSymbol.reductionType();
			this.conductorExponent = conductorExponent;
			this.localIndex = localIndex;
			init();
		}

		private void init() {
			reductionMap = new MathMap<>() {
				@Override
				public ProjectivePoint<S> evaluate(ProjectivePoint<T> t) {
					if (t.equals(translation.getDomain().neutral())) {
						return new ProjectivePoint<>(field.residueField(), field.residueField().zero(),
								field.residueField().one(), field.residueField().zero());
					}
					ProjectivePoint<T> translated = translation.evaluate(t);
					return new ProjectivePoint<>(field.residueField(),
							field.reduceInteger(translated.getDehomogenisedCoord(1, 3)),
							field.reduceInteger(translated.getDehomogenisedCoord(2, 3)), field.residueField().one());
				}
			};
		}

		public String getKodairaSymbolString() {
			switch (kodairaSymbol) {
			case I0:
			case I1:
			case I2:
			case II:
			case III:
			case IV:
				return kodairaSymbol.name();
			case Iv:
				return "I" + getKodairaSymbolValue();
			case mIv:
				return getKodairaSymbolMultiplicity() + "I" + getKodairaSymbolValue();
			case I0star:
				return "I0*";
			case Ivstar:
				return "I" + getKodairaSymbolValue() + "*";
			case IIstar:
				return "II*";
			case IIIstar:
				return "III*";
			case IVstar:
				return "IV*";
			default:
				throw new ArithmeticException("Unknown Kodaira Symbol!");
			}
		}

		public KodairaSymbol getKodairaSymbol() {
			return kodairaSymbol;
		}

		public int getKodairaSymbolValue() {
			return value;
		}

		public int getKodairaSymbolMultiplicity() {
			return multiplicity;
		}

		public ReductionType getReductionType() {
			return reductionType;
		}

		public boolean hasGoodReduction() {
			return reductionType == ReductionType.GOOD_REDUCTION;
		}

		public int getConductorExponent() {
			return conductorExponent;
		}

		public int getLocalIndex() {
			return localIndex;
		}

		public GenericProjectiveScheme<S> getReducedScheme() {
			if (reducedScheme == null) {
				if (kodairaSymbol.equals(KodairaSymbol.I0)) {
					reducedScheme = reducedCurve.asGenericProjectiveScheme();
				} else {
					PolynomialRing<S> projectivePolynomialRing = AbstractPolynomialRing
							.getPolynomialRing(field.residueField(), 3, Monomial.GREVLEX);
					Map<Monomial, S> coeffs = new TreeMap<>();
					coeffs.put(projectivePolynomialRing.getMonomial(new int[] { 3, 0, 0 }), field.residueField().one());
					coeffs.put(projectivePolynomialRing.getMonomial(new int[] { 0, 2, 1 }),
							field.residueField().getInteger(-1));
					coeffs.put(projectivePolynomialRing.getMonomial(new int[] { 1, 1, 1 }),
							field.reduceInteger(field.negative(translation.getRange().getA1())));
					coeffs.put(projectivePolynomialRing.getMonomial(new int[] { 2, 0, 1 }),
							field.reduceInteger(translation.getRange().getA2()));
					coeffs.put(projectivePolynomialRing.getMonomial(new int[] { 0, 1, 2 }),
							field.reduceInteger(field.negative(translation.getRange().getA3())));
					coeffs.put(projectivePolynomialRing.getMonomial(new int[] { 1, 0, 2 }),
							field.reduceInteger(translation.getRange().getA4()));
					coeffs.put(projectivePolynomialRing.getMonomial(new int[] { 0, 0, 3 }),
							field.reduceInteger(translation.getRange().getA6()));
					reducedScheme = new GenericProjectiveScheme<>(field.residueField(), projectivePolynomialRing,
							Collections.singletonList(projectivePolynomialRing.getPolynomial(coeffs)));
				}
			}
			return reducedScheme;
		}

		public EllipticCurve<T> getCurve() {
			return translation.getDomain();
		}

		public EllipticCurve<T> getMinimalModel() {
			return translation.getRange();
		}

		public Isomorphism<T> getTranslation() {
			return translation;
		}

		public EllipticCurve<S> getReducedCurve() {
			return reducedCurve;
		}

		public MathMap<ProjectivePoint<T>, ProjectivePoint<S>> getReductionMap() {
			return reductionMap;
		}
	}

	public <S extends Element<S>> TatesAlgorithmResult<T, S> tatesAlgorithm(DiscreteValuationField<T, S> field) {
		if (!field.isInteger(a1) || !field.isInteger(a2) || !field.isInteger(a3) || !field.isInteger(a4)
				|| !field.isInteger(a6)) {
			throw new ArithmeticException("Not an integer curve!");
		}
		Isomorphism<T> translation = this.identity();
		EllipticCurve<T> model = translation.getRange();
		Field<S> reduction = field.residueField();
		while (true) {
			if (field.valuation(model.discriminant()).equals(Value.ZERO)) {
				return new TatesAlgorithmResult<>(field, translation, KodairaSymbol.I0, 0, 1);
			}
			T xShift;
			T yShift;
			if (reduction.characteristic().equals(BigInteger.TWO)) {
				if (field.valuation(model.getB2()).compareTo(Value.ONE) >= 0) {
					xShift = field.liftToInteger(reduction.characteristicRoot(field.reduceInteger(model.getA4())));
					S reducedXShift = field.reduceInteger(xShift);
					yShift = field.liftToInteger(reduction.characteristicRoot(
							reduction.add(reduction.multiply(reducedXShift, field.reduceInteger(model.getA4())),
									reduction.multiply(reducedXShift, reducedXShift,
											field.reduceInteger(model.getA2())),
									reduction.power(reducedXShift, 3), field.reduceInteger(model.getA6()))));
				} else {
					xShift = field.liftToInteger(
							reduction.divide(field.reduceInteger(model.getA3()), field.reduceInteger(model.getA1())));
					yShift = field.liftToInteger(reduction.divide(
							reduction.add(field.reduceInteger(model.getA4()),
									reduction.power(field.reduceInteger(xShift), 2)),
							field.reduceInteger(model.getA1())));
				}
			} else if (reduction.characteristic().equals(BigInteger.valueOf(3))) {
				if (field.valuation(model.getB2()).compareTo(Value.ONE) >= 0) {
					xShift = field.liftToInteger(
							reduction.negative(reduction.characteristicRoot(field.reduceInteger(model.getB6()))));
				} else {
					xShift = field
							.liftToInteger(reduction.divide(reduction.negative(field.reduceInteger(model.getB4())),
									field.reduceInteger(model.getB2())));
				}
				yShift = field.round(field.add(field.multiply(xShift, model.getA1()), model.getA3()), 1);
			} else {
				if (field.valuation(model.getC4()).compareTo(Value.ONE) >= 0) {
					xShift = field.multiply(-1, field.liftToInteger(reduction.inverse(reduction.getInteger(12))),
							model.getB2());
				} else {
					xShift = field.multiply(-1,
							field.liftToInteger(
									reduction.inverse(reduction.multiply(12, field.reduceInteger(model.getC4())))),
							field.add(model.getC6(), field.multiply(model.getB2(), model.getC4())));
				}
				yShift = field.multiply(-1, field.liftToInteger(reduction.inverse(reduction.getInteger(2))),
						field.add(field.multiply(xShift, model.getA1()), model.getA3()));
				xShift = field.round(xShift, 1);
				yShift = field.round(yShift, 1);
			}
			Isomorphism<T> shift = new Isomorphism<T>(model, field.one(), xShift, field.zero(), yShift);
			translation = concatIsomorphisms(translation, shift);
			model = translation.getRange();
			if (field.valuation(model.getA3()).compareTo(Value.ZERO) <= 0
					|| field.valuation(model.getA4()).compareTo(Value.ZERO) <= 0
					|| field.valuation(model.getA6()).compareTo(Value.ZERO) <= 0) {
				throw new ArithmeticException("Coordinate change not successful!");
			}
			if (field.valuation(model.getC4()).compareTo(Value.ONE) < 0) {
				KodairaSymbol symbol = KodairaSymbol.Iv;
				int conductorExponent = 1;
				List<S> indexCoeffs = new ArrayList<>();
				indexCoeffs.add(field.reduceInteger(field.negative(model.getA2())));
				indexCoeffs.add(field.reduceInteger(model.getA1()));
				indexCoeffs.add(reduction.one());
				UnivariatePolynomial<S> indexPolynomial = reduction.getUnivariatePolynomialRing()
						.getPolynomial(indexCoeffs);
				int localIndex;
				int discriminantValue = field.valuation(model.discriminant()).value();
				if (reduction.hasRoots(indexPolynomial)) {
					localIndex = discriminantValue;
				} else if (discriminantValue % 2 == 0) {
					localIndex = 2;
				} else {
					localIndex = 1;
				}
				return new TatesAlgorithmResult<>(field, translation, symbol, discriminantValue, conductorExponent,
						localIndex);
			}
			if (field.valuation(model.getA6()).compareTo(new Value(2)) < 0) {
				KodairaSymbol symbol = KodairaSymbol.II;
				int conductorExponent = field.valuation(model.discriminant()).value();
				int localIndex = 1;
				return new TatesAlgorithmResult<>(field, translation, symbol, conductorExponent, localIndex);
			}
			if (field.valuation(model.getB8()).compareTo(new Value(3)) < 0) {
				KodairaSymbol symbol = KodairaSymbol.III;
				int conductorExponent = field.valuation(model.discriminant()).value() - 1;
				int localIndex = 2;
				return new TatesAlgorithmResult<>(field, translation, symbol, conductorExponent, localIndex);
			}
			if (field.valuation(model.getB6()).compareTo(new Value(3)) < 0) {
				KodairaSymbol symbol = KodairaSymbol.IV;
				int conductorExponent = field.valuation(model.discriminant()).value() - 2;
				List<S> indexCoeffs = new ArrayList<>();
				indexCoeffs.add(field.reduceInteger(
						field.negative(field.divide(model.getA6(), field.power(field.uniformizer(), 2)))));
				indexCoeffs.add(field.reduceInteger(field.divide(model.getA3(), field.uniformizer())));
				indexCoeffs.add(reduction.one());
				UnivariatePolynomial<S> indexPolynomial = reduction.getUnivariatePolynomialRing()
						.getPolynomial(indexCoeffs);
				int localIndex = reduction.hasRoots(indexPolynomial) ? 3 : 1;
				return new TatesAlgorithmResult<>(field, translation, symbol, conductorExponent, localIndex);
			}
			UnivariatePolynomialRing<S> r = reduction.getUnivariatePolynomialRing();
			T xyShift;
			yShift = null;
			if (reduction.characteristic().equals(BigInteger.TWO)) {
				xyShift = field.liftToInteger(reduction.characteristicRoot(field.reduceInteger(model.getA2())));
				yShift = field.multiply(field.uniformizer(), field.liftToInteger(reduction.characteristicRoot(
						field.reduceInteger(field.divide(model.getA6(), field.power(field.uniformizer(), 2))))));
			} else {
				xyShift = field.round(field.multiply(-1, model.getA1(),
						field.liftToInteger(reduction.inverse(reduction.getInteger(2)))), 1);
				yShift = field.multiply(-1, model.getA3(),
						field.liftToInteger(reduction.inverse(reduction.getInteger(2))));
			}
			shift = new Isomorphism<T>(model, field.one(), field.zero(), xyShift, yShift);
			translation = concatIsomorphisms(translation, shift);
			model = translation.getRange();
			if (field.valuation(model.getA1()).compareTo(Value.ZERO) <= 0
					|| field.valuation(model.getA2()).compareTo(Value.ZERO) <= 0
					|| field.valuation(model.getA3()).compareTo(Value.ONE) <= 0
					|| field.valuation(model.getA4()).compareTo(Value.ONE) <= 0
					|| field.valuation(model.getA6()).compareTo(new Value(2)) <= 0) {
				throw new ArithmeticException("Change of variables failed!");
			}
			S b = field.reduceInteger(field.divide(model.getA2(), field.uniformizer()));
			S c = field.reduceInteger(field.divide(model.getA4(), field.power(field.uniformizer(), 2)));
			S d = field.reduceInteger(field.divide(model.getA6(), field.power(field.uniformizer(), 3)));
			S x = reduction.subtract(reduction.multiply(3, c), reduction.multiply(b, b));
			List<S> coeffsRhs = new ArrayList<>();
			coeffsRhs.add(d);
			coeffsRhs.add(c);
			coeffsRhs.add(b);
			coeffsRhs.add(reduction.one());
			UnivariatePolynomial<S> reducedRhs = r.getPolynomial(coeffsRhs);
			FactorizationResult<Polynomial<S>, S> reducedRhsFactors = reduction.factorization(reducedRhs);
			if (reducedRhsFactors.squareFree()) {
				KodairaSymbol symbol = KodairaSymbol.I0star;
				int conductorExponent = field.valuation(model.discriminant()).value() - 4;
				int localIndex = 1 + reduction.roots(reducedRhs).size();
				return new TatesAlgorithmResult<>(field, translation, symbol, conductorExponent, localIndex);
			}
			if (!x.equals(reduction.zero())) {
				KodairaSymbol symbol = KodairaSymbol.Ivstar;
				xShift = null;
				if (reduction.characteristic().equals(BigInteger.TWO)) {
					xShift = field.liftToInteger(reduction.characteristicRoot(c));
				} else if (reduction.characteristic().equals(BigInteger.valueOf(3))) {
					xShift = field.liftToInteger(reduction.multiply(b, c));
				} else {
					xShift = field.liftToInteger(
							reduction.divide(reduction.subtract(reduction.multiply(b, c), reduction.multiply(9, d)),
									reduction.multiply(2, x)));
				}
				xShift = field.multiply(field.uniformizer(), xShift);
				shift = new Isomorphism<T>(model, field.one(), xShift, field.zero(), field.zero());
				translation = concatIsomorphisms(translation, shift);
				model = translation.getRange();
				int value = 1;
				T mx = field.power(field.uniformizer(), 2);
				T my = field.power(field.uniformizer(), 2);
				int localIndex = 0;
				while (localIndex == 0) {
					T xa2 = field.divide(model.getA2(), field.uniformizer());
					T xa3 = field.divide(model.getA3(), my);
					T xa4 = field.divide(model.getA4(), field.multiply(field.uniformizer(), mx));
					T xa6 = field.divide(model.getA6(), field.multiply(mx, my));
					if (field.valuation(xa2).compareTo(Value.ZERO) < 0 || field.valuation(xa3).compareTo(Value.ZERO) < 0
							|| field.valuation(xa4).compareTo(Value.ZERO) < 0
							|| field.valuation(xa6).compareTo(Value.ZERO) < 0) {
						throw new ArithmeticException("Coordinate change not successful");
					}
					if (field.valuation(field.add(field.multiply(xa3, xa3), field.multiply(4, xa6)))
							.compareTo(Value.ONE) < 0) {
						List<S> indexCoeffs = new ArrayList<>();
						indexCoeffs.add(field.reduceInteger(field.negative(xa6)));
						indexCoeffs.add(field.reduceInteger(xa3));
						indexCoeffs.add(reduction.one());
						UnivariatePolynomial<S> indexPolynomial = reduction.getUnivariatePolynomialRing()
								.getPolynomial(indexCoeffs);
						localIndex = reduction.hasRoots(indexPolynomial) ? 4 : 2;
					} else {
						yShift = null;
						if (reduction.characteristic().equals(BigInteger.TWO)) {
							yShift = field.multiply(my,
									field.liftToInteger(reduction.characteristicRoot(field.reduceInteger(xa6))));
						} else {
							yShift = field.multiply(my, field.liftToInteger(reduction
									.negative(reduction.divide(field.reduceInteger(xa3), reduction.getInteger(2)))));
						}
						translation = concatIsomorphisms(translation, new Isomorphism<T>(translation.getRange(),
								field.one(), field.zero(), field.zero(), yShift));
						model = translation.getRange();
						my = field.multiply(field.uniformizer(), my);
						value++;
						xa2 = field.divide(model.getA2(), field.uniformizer());
						xa3 = field.divide(model.getA3(), my);
						xa4 = field.divide(model.getA4(), field.multiply(field.uniformizer(), mx));
						xa6 = field.divide(model.getA6(), field.multiply(mx, my));
						if (field.valuation(xa2).compareTo(Value.ZERO) < 0
								|| field.valuation(xa3).compareTo(Value.ZERO) < 0
								|| field.valuation(xa4).compareTo(Value.ZERO) < 0
								|| field.valuation(xa6).compareTo(Value.ZERO) < 0) {
							throw new ArithmeticException("Coordinate change not successful");
						}
						if (field.valuation(field.subtract(field.multiply(xa4, xa4), field.multiply(4, xa2, xa6)))
								.compareTo(Value.ONE) < 0) {
							List<S> indexCoeffs = new ArrayList<>();
							indexCoeffs.add(field.reduceInteger(xa6));
							indexCoeffs.add(field.reduceInteger(xa4));
							indexCoeffs.add(field.reduceInteger(xa2));
							UnivariatePolynomial<S> indexPolynomial = reduction.getUnivariatePolynomialRing()
									.getPolynomial(indexCoeffs);
							localIndex = reduction.hasRoots(indexPolynomial) ? 4 : 2;
						} else {
							xShift = null;
							if (reduction.characteristic().equals(BigInteger.TWO)) {
								if (field.valuation(xa2).equals(Value.ZERO)) {
									xShift = field.multiply(mx, field.liftToInteger(reduction.characteristicRoot(
											reduction.divide(field.reduceInteger(xa6), field.reduceInteger(xa2)))));
								} else {
									xShift = field.zero();
								}
							} else {
								xShift = field.multiply(mx,
										field.liftToInteger(
												reduction.divide(reduction.negative(field.reduceInteger(xa4)),
														reduction.multiply(2, field.reduceInteger(xa2)))));
							}
							translation = concatIsomorphisms(translation, new Isomorphism<T>(translation.getRange(),
									field.one(), xShift, field.zero(), field.zero()));
							model = translation.getRange();
							mx = field.multiply(field.uniformizer(), mx);
							value++;
						}
					}
				}
				int conductorExponent = field.valuation(model.discriminant()).value() - 4 - value;
				return new TatesAlgorithmResult<>(field, translation, symbol, value, conductorExponent, localIndex);
			}
			xShift = null;
			if (reduction.characteristic().equals(BigInteger.valueOf(3))) {
				xShift = field.negative(field.liftToInteger(reduction.characteristicRoot(d)));
			} else {
				xShift = field.liftToInteger(reduction.divide(reduction.negative(b), reduction.getInteger(3)));
			}
			xShift = field.multiply(field.uniformizer(), xShift);
			shift = new Isomorphism<T>(model, field.one(), xShift, field.zero(), field.zero());
			translation = concatIsomorphisms(translation, shift);
			model = translation.getRange();
			S x3 = field.reduceInteger(field.divide(model.getA3(), field.power(field.uniformizer(), 2)));
			S x6 = field.reduceInteger(field.divide(model.getA6(), field.power(field.uniformizer(), 4)));
			if (!reduction.add(reduction.multiply(x3, x3), reduction.multiply(4, x6)).equals(reduction.zero())) {
				KodairaSymbol symbol = KodairaSymbol.IVstar;
				int conductorExponent = field.valuation(model.discriminant()).value() - 6;
				List<S> indexCoeffs = new ArrayList<>();
				indexCoeffs.add(reduction.negative(x6));
				indexCoeffs.add(x3);
				indexCoeffs.add(reduction.one());
				UnivariatePolynomial<S> indexPolynomial = reduction.getUnivariatePolynomialRing()
						.getPolynomial(indexCoeffs);
				int localIndex = reduction.hasRoots(indexPolynomial) ? 3 : 1;
				return new TatesAlgorithmResult<>(field, translation, symbol, conductorExponent, localIndex);
			}
			yShift = null;
			if (reduction.characteristic().equals(BigInteger.TWO)) {
				yShift = field.liftToInteger(reduction.characteristicRoot(x6));
			} else {
				yShift = field.liftToInteger(reduction.divide(x3, reduction.getInteger(2)));
			}
			yShift = field.multiply(-1, field.power(field.uniformizer(), 2), yShift);
			shift = new Isomorphism<T>(model, field.one(), field.zero(), field.zero(), yShift);
			translation = concatIsomorphisms(translation, shift);
			model = translation.getRange();
			if (field.valuation(model.getA4()).compareTo(new Value(4)) < 0) {
				KodairaSymbol symbol = KodairaSymbol.IIIstar;
				int conductorExponent = field.valuation(model.discriminant()).value() - 7;
				int localIndex = 2;
				return new TatesAlgorithmResult<>(field, translation, symbol, conductorExponent, localIndex);

			}
			if (field.valuation(model.getA6()).compareTo(new Value(6)) < 0) {
				KodairaSymbol symbol = KodairaSymbol.IIstar;
				int conductorExponent = field.valuation(model.discriminant()).value() - 8;
				int localIndex = 1;
				return new TatesAlgorithmResult<>(field, translation, symbol, conductorExponent, localIndex);
			}
			shift = new Isomorphism<T>(model, field.uniformizer(), field.zero(), field.zero(), field.zero());
			translation = concatIsomorphisms(translation, shift);
			model = translation.getRange();
		}
	}

	public static class StableModel<T extends Element<T>, I extends Element<I>, S extends Element<S>, B extends Element<B>, IB extends Element<IB>, SB extends Element<SB>, E extends AlgebraicExtensionElement<B, E>, IE extends Element<IE>, R extends Element<R>, RE extends AlgebraicExtensionElement<R, RE>, RFE extends FieldExtension<R, RE, RFE>, DFE extends DiscreteValuationFieldExtension<B, SB, E, DFE, R, RE, RFE>, FE extends GlobalFieldExtension<B, IB, SB, E, IE, R, RE, RFE, DFE, FE>> {
		private GlobalField<T, I, S> field;
		private ExtensionOfGlobalField<T, I, S, B, IB, SB, E, IE, R, RE, RFE, DFE, FE> fieldExtension;
		private EllipticCurve<T> curve;
		private EllipticCurve<E> stableModel;
	}

	public <I extends Element<I>, S extends Element<S>> StableModel<T, I, S, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?> stableModel(
			GlobalField<T, I, S> field) {
		return stableModel(field, field.getGlobalFieldExtension(field.getUnivariatePolynomialRing().getVar()));
	}

	public <I extends Element<I>, S extends Element<S>, B extends Element<B>, IB extends Element<IB>, SB extends Element<SB>, E extends AlgebraicExtensionElement<B, E>, IE extends Element<IE>, R extends Element<R>, RE extends AlgebraicExtensionElement<R, RE>, RFE extends FieldExtension<R, RE, RFE>, DFE extends DiscreteValuationFieldExtension<B, SB, E, DFE, R, RE, RFE>, FE extends GlobalFieldExtension<B, IB, SB, E, IE, R, RE, RFE, DFE, FE>> StableModel<T, I, S, B, IB, SB, E, IE, R, RE, RFE, DFE, FE> stableModel(
			GlobalField<T, I, S> field,
			ExtensionOfGlobalField<T, I, S, B, IB, SB, E, IE, R, RE, RFE, DFE, FE> trivialExtension) {
		Integers z = Integers.z();
		FE extension = trivialExtension.getExtension();
		UnivariatePolynomialRing<E> polynomials = extension.getUnivariatePolynomialRing();
		EllipticCurve<E> extendedCurve = extendBaseField(trivialExtension.asExtension()).getCurve();
		extendedCurve = extendedCurve.minimalModel(extension).get().getRange();
		DedekindRingExtension<B, IB, SB, E, IE, R, RE, RFE, DFE, FE> integers = extension.ringOfIntegers();
		Ideal<IE> discriminantIdeal = integers
				.getIdeal(Collections.singletonList(integers.asInteger(extendedCurve.discriminant())));
		FactorizationResult<Ideal<IE>, Ideal<IE>> discriminantFactors = integers.idealFactorization(discriminantIdeal);
		FieldEmbedding<B, E, FE> embedding = extension.getEmbeddedExtension(polynomials.getVar());
		for (Ideal<IE> prime : discriminantFactors.primeFactors()) {
			DFE localField = integers.localize(prime).localField();
			UnivariatePolynomial<E> minimalPolynomial;
			if (!localField.residueField().characteristic().equals(BigInteger.TWO)
					&& !localField.residueField().characteristic().equals(BigInteger.valueOf(3))) {
				IntE value = z.getInteger(localField.valuation(extendedCurve.discriminant()).value());
				IntE gcd = z.gcd(value, z.getInteger(12));
				int power = 12 / gcd.intValueExact();
				if (power == 1) {
					continue;
				}
				minimalPolynomial = polynomials.toUnivariate(polynomials.subtract(polynomials.getVarPower(power),
						polynomials.getEmbedding(localField.uniformizer())));
			} else if (localField.residueField().characteristic().equals(BigInteger.TWO)) {
				if (localField.valuation(extendedCurve.jInvariant()).compareTo(Value.ZERO) > 0) {

				} else {
					E denominator = extension
							.inverse(extension.subtract(extendedCurve.jInvariant(), extension.getInteger(1728)));
					EllipticCurve<E> otherModel = new EllipticCurve<>(extension, extension.one(), extension.zero(),
							extension.zero(), extension.multiply(-36, denominator), extension.negative(denominator));
				}
			} else {

			}
			minimalPolynomial = null;
			embedding = new FieldEmbedding<>(embedding, extension.getEmbeddedExtension(minimalPolynomial));
			extension = embedding.getField();
			polynomials = extension.getUnivariatePolynomialRing();
			extendedCurve = extendBaseField(trivialExtension.extendFurther(embedding).asExtension()).getCurve();
			integers = extension.ringOfIntegers();
		}
		extendedCurve = extendedCurve.minimalModel(extension).get().getRange();
		EllipticCurve<E> overFieldExtension = extendBaseField(trivialExtension.asExtension()).getCurve();
		return null;
	}

	public ProjectiveMorphism<T> xCover() {
		if (xCover == null) {
			List<Polynomial<T>> cover = new ArrayList<>();
			cover.add(projectiveRing.getVar(1));
			cover.add(projectiveRing.getVar(3));
			xCover = new ProjectiveMorphism<>(asGenericProjectiveScheme(),
					new ProjectiveLine<>(field).asGenericProjectiveScheme(), cover);
		}
		return xCover;
	}

	public ProjectivePoint<T> kummerCoordinate(ProjectivePoint<T> point) {
		return xCover().evaluate(point);
	}

	// file:///Users/sophie/Downloads/Computing_Canonical_Heights_on_Elliptic_Curves_in_.pdf
	private Real nonArchimedeanLocalHeight(EmbeddedNumberField<Ext, CompletedNumberField> field,
			ProjectivePoint<T> point) {
		Reals r = field.getReals();
		Rationals q = Rationals.q();
		Integers z = Integers.z();
		if (point.equals(neutral())) {
			return r.zero();
		}
		NFE x = (NFE) point.getDehomogenisedCoord(1, 3);
		Real lambda = r.log(r.max(field.value(x), r.one()));
		int bound = field.embeddingField().valuation(field.embedding((NFE) discriminant())).value();
		if (bound <= 1) {
			return lambda;
		}
		Ext b2 = field.embedding((NFE) getB2());
		Ext b4 = field.embedding((NFE) getB4());
		Ext b6 = field.embedding((NFE) getB6());
		Ext b8 = field.embedding((NFE) getB8());
		PolynomialRing<Ext> ring = AbstractPolynomialRing.getPolynomialRing(field.embeddingField(), Monomial.GREVLEX,
				new String[] { "W", "Z" });
		List<Polynomial<Ext>> deltaList = new ArrayList<>();
		deltaList.add(ring.add(ring.getVarPower(1, 4),
				ring.multiply(field.embeddingField().multiply(-1, b4), ring.getVarPower(1, 2), ring.getVarPower(2, 2)),
				ring.multiply(field.embeddingField().multiply(-2, b6), ring.getVar(1), ring.getVarPower(2, 3)),
				ring.multiply(field.embeddingField().multiply(-1, b8), ring.getVarPower(2, 4))));
		deltaList.add(ring.add(ring.multiply(4, ring.getVarPower(1, 3), ring.getVar(2)),
				ring.multiply(b2, ring.getVarPower(1, 2), ring.getVarPower(2, 2)),
				ring.multiply(field.embeddingField().multiply(2, b4), ring.getVar(1), ring.getVarPower(2, 3)),
				ring.multiply(b6, ring.getVarPower(2, 4))));
		Vector<Polynomial<Ext>> delta = new Vector<>(deltaList);
		int m = r.divide(r.log(r.divide(r.power(r.getInteger(bound), 3), r.getInteger(3))), r.log(r.getInteger(4)))
				.roundDown().intValueExact();
		int accuracy = (m + 1) * bound + 1;
		accuracy = Math.max(accuracy, 3);
		accuracy *= 2;
		if (field.embeddingField().getAccuracy() < accuracy) {
			List<EmbeddedNumberField<Ext, CompletedNumberField>> embeddings = field.numberField()
					.padicEmbeddings(field.embeddingField().getBaseField().withAccuracy(accuracy));
			for (EmbeddedNumberField<Ext, CompletedNumberField> embeddedField : embeddings) {
				if (embeddedField.embeddingField().exactIdeal().equals(field.embeddingField().exactIdeal())) {
					field = embeddedField;
					break;
				}
			}
		}
		Fraction mu0 = q.zero();
		FreeModule<Ext> module = new FreeModule<>(field.embeddingField(), 2);
		Vector<Ext> kummer = new Vector<>(
				field.embeddingField().round(field.embedding((NFE) point.getDehomogenisedCoord(1, 3)), accuracy),
				field.embeddingField().one());
		for (int n = 0; n <= m; n++) {
			kummer = field.embeddingField().ringOfIntegers().roundVector(ring.evaluate(delta, kummer), accuracy);
			int minValue = field.embeddingField().valuation(kummer.get(1))
					.min(field.embeddingField().valuation(kummer.get(2))).value();
			if (minValue == 0) {
				return r.subtract(lambda,
						r.multiply(r.log(r.getInteger(field.embeddingField().residueField().getNumberOfElements())),
								r.getFraction(q.divide(mu0, q.getInteger(field.embeddingField().degree())))));
			}
			mu0 = q.add(mu0, q.getFraction(z.getInteger(minValue), z.power(z.getInteger(4), n + 1)));
			kummer = module.scalarMultiply(
					field.embeddingField().power(field.embeddingField().uniformizer(), -minValue), kummer);
		}
		Fraction lowerBound = mu0;
		Fraction upperBound = q.add(mu0, q.inverse(q.getInteger(bound * bound)));
		Fraction middle = q.divide(q.add(lowerBound, upperBound), q.getInteger(2));
		Iterator<Fraction> muIterator = MiscAlgorithms.continuedFractionApproximation(q, middle, new MathMap<>() {
			@Override
			public IntE evaluate(Fraction t) {
				return t.roundDown();
			}
		});
		Fraction mu;
		do {
			mu = muIterator.next();
		} while (lowerBound.compareTo(mu) > 0 || upperBound.compareTo(mu) < 0);
		if (mu.getDenominator().compareTo(z.getInteger(bound)) > 0) {
			throw new ArithmeticException("Did not find mu!");
		}
		return r.subtract(lambda,
				r.multiply(r.log(r.getInteger(field.embeddingField().residueField().getNumberOfElements())),
						r.getFraction(q.divide(mu0, q.getInteger(field.embeddingField().degree())))));
	}

	// https://www.ams.org/journals/mcom/1988-51-183/S0025-5718-1988-0942161-4/S0025-5718-1988-0942161-4.pdf
	@SuppressWarnings("unchecked")
	private <S extends Element<S>, F extends ValueField<S>> Real archimedeanLocalHeight(EmbeddedNumberField<S, F> field,
			ProjectivePoint<T> point) {
		Reals r = field.getReals();
		if (point.equals(neutral())) {
			return r.zero();
		}
		S[] b2 = (S[]) Array.newInstance(field.embeddingField().zero().getClass(), 2);
		S[] b4 = (S[]) Array.newInstance(field.embeddingField().zero().getClass(), 2);
		S[] b6 = (S[]) Array.newInstance(field.embeddingField().zero().getClass(), 2);
		S[] b8 = (S[]) Array.newInstance(field.embeddingField().zero().getClass(), 2);
		b2[0] = field.embedding((NFE) getB2());
		b4[0] = field.embedding((NFE) getB4());
		b6[0] = field.embedding((NFE) getB6());
		b8[0] = field.embedding((NFE) getB8());
		S x = field.embedding((NFE) point.getDehomogenisedCoord(1, 3));
		b2[1] = field.embeddingField().add(b2[0], field.embeddingField().getInteger(-12));
		b4[1] = field.embeddingField().add(b4[0], field.embeddingField().multiply(-1, b2[0]),
				field.embeddingField().getInteger(6));
		b6[1] = field.embeddingField().add(b6[0], field.embeddingField().multiply(-2, b4[0]), b2[0],
				field.embeddingField().getInteger(-4));
		b8[1] = field.embeddingField().add(b8[0], field.embeddingField().multiply(-3, b6[0]),
				field.embeddingField().multiply(3, b4[0]), field.embeddingField()
						.add(field.embeddingField().multiply(-1, b2[0]), field.embeddingField().getInteger(3)));
		S t;
		int beta;
		if (field.embeddingField().value(x).compareTo(r.getDouble(0.5)) >= 0) {
			t = field.embeddingField().inverse(x);
			beta = 0;
		} else {
			t = field.embeddingField().inverse(field.embeddingField().add(x, field.embeddingField().one()));
			beta = 1;
		}
		int n = 0;
		Real lambda = r.multiply(r.getDouble(-0.5), r.log(field.embeddingField().value(t)));
		Real mu = r.zero();
		Real prevMu;
		do {
			prevMu = mu;
			S w = field.embeddingField().add(field.embeddingField().multiply(4, t),
					field.embeddingField().multiply(b2[beta], t, t),
					field.embeddingField().multiply(2, b4[beta], field.embeddingField().power(t, 3)),
					field.embeddingField().multiply(b6[beta], field.embeddingField().power(t, 4)));
			S z = field.embeddingField().add(field.embeddingField().one(),
					field.embeddingField().multiply(-1, b4[beta], t, t),
					field.embeddingField().multiply(-2, b6[beta], field.embeddingField().power(t, 3)),
					field.embeddingField().multiply(-1, b8[beta], field.embeddingField().power(t, 4)));
			if (field.embeddingField().value(w).compareTo(r.multiply(2, field.embeddingField().value(z))) <= 0) {
				mu = r.add(r.multiply(r.getPowerOfTwo(-2 * n), r.log(field.embeddingField().value(z))), mu);
				t = field.embeddingField().divide(w, z);
			} else {
				mu = r.add(
						r.multiply(
								r.getPowerOfTwo(-2
										* n),
								r.log(field.embeddingField()
										.value(field.embeddingField().add(z,
												field.embeddingField().multiply(field.embeddingField()
														.power(field.embeddingField().getInteger(-1), beta), w))))),
						mu);
				t = field.embeddingField().divide(w, field.embeddingField().add(z, field.embeddingField()
						.multiply(field.embeddingField().power(field.embeddingField().getInteger(-1), beta), w)));
				beta = 1 - beta;
			}
			n++;
		} while (!r.close(mu, prevMu));
		return /* r.add( */r.multiply(2,
				r.add(lambda, r.multiply(r.getPowerOfTwo(-3), mu)));/*
																	 * , r.multiply(r.inverse(r.getInteger(6)),r.log(
																	 * field.value((NFE) discriminant()))));
																	 */
	}

	@SuppressWarnings("unchecked")
	public <S extends Element<S>, F extends ValueField<S>> Real localHeight(EmbeddedNumberField<S, F> field,
			ProjectivePoint<T> point) {
		if (!(this.field instanceof NumberField)) {
			throw new ArithmeticException("Only defined over number field!");
		}
		if (field.embeddingField() instanceof CompletedNumberField) {
			return nonArchimedeanLocalHeight((EmbeddedNumberField<Ext, CompletedNumberField>) field, point);
		}
		return archimedeanLocalHeight(field, point);
	}

	public Real height(ProjectivePoint<T> point) {
		if (!(field instanceof NumberField)) {
			throw new ArithmeticException("Only defined over number field!");
		}
		NumberField nf = (NumberField) field;
		Reals r = nf.minkowskiEmbeddingSpace().getValueField();
		if (point.equals(neutral())) {
			return r.zero();
		}
		NFE discriminant = (NFE) discriminant();
		Set<Ideal<NFE>> primes = new TreeSet<>();
		primes.addAll(nf.maximalOrder().idealFactorization(nf.maximalOrder().getIdeal(discriminant)).primeFactors());
		NFE x = (NFE) point.getDehomogenisedCoord(1, 3);
		primes.addAll(nf.maximalOrder()
				.idealFactorization(nf.maximalOrder().getIdeal(nf.maximalOrder().getDenominator(x))).primeFactors());
		Map<IntE, Integer> intPrimes = new TreeMap<>();
		Real result = r.zero();
		for (Ideal<NFE> prime : primes) {
			LocalizedNumberField local = nf.maximalOrder().localizeAndQuotient(prime);
			int bound = Math.max(local.valuation((NFE) discriminant()).value(), 1);
			int m = r.divide(r.log(r.divide(r.power(r.getInteger(bound), 3), r.getInteger(3))), r.log(r.getInteger(4)))
					.roundDown().intValueExact();
			int accuracy = (m + 1) * bound + 1;
			accuracy = Math.max(accuracy, 3);
			accuracy *= 3;
			IntE intPrime = ((NumberFieldIdeal) prime).intGenerator();
			if (!intPrimes.containsKey(intPrime) || intPrimes.get(intPrime) < accuracy) {
				intPrimes.put(intPrime, accuracy);
			}
		}
		for (IntE prime : intPrimes.keySet()) {
			PAdicField padic = new PAdicField(prime, intPrimes.get(prime));
			for (EmbeddedNumberField<Ext, CompletedNumberField> local : nf.padicEmbeddings(padic)) {
				result = r.add(r.multiply(local.embeddingField().degree(), localHeight(local, point)), result);
			}
		}
		for (EmbeddedNumberField<Real, Reals> local : nf.realEmbeddings()) {
			result = r.add(localHeight(local, point), result);
		}
		for (EmbeddedNumberField<ComplexNumber, Complex> local : nf.complexEmbeddings()) {
			result = r.add(r.multiply(2, localHeight(local, point)), result);
		}
		return r.divide(result, r.getInteger(nf.degree()));
	}

	public Real neronTatePairing(ProjectivePoint<T> p, ProjectivePoint<T> q) {
		Reals r = ((NumberField) field).minkowskiEmbeddingSpace().getValueField();
		return r.divide(r.subtract(height(add(p, q)), r.add(height(p), height(q))), r.getInteger(2));
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
			if (!field.isFinite()) {
				throw new ArithmeticException("Field is not finite!");
			}
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
				} while (tries < 32 && multiply(subPrimePower, torsionPoint1).equals(neutral()));
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
				if (z.isDivisible(order, z.power(prime, 2 * multiplicity + multiplicationOffset))
						&& z.isDivisible(z.getInteger(field.getNumberOfUnits()), torsion)) {
					ProjectivePoint<T> torsionPoint2;
					tries = 0;
					T weilPairing = field.one();
					do {
						torsionPoint2 = multiply(multiplier, getRandomElement());
						tries++;
						if (multiply(subPrimePower, torsionPoint2).equals(neutral())) {
							continue;
						}
						while (!multiply(torsion, torsionPoint2).equals(neutral())) {
							torsionPoint2 = multiply(prime, torsionPoint2);
						}
						weilPairing = weilPairing(torsion.getValue(), torsionPoint1, torsionPoint2);
					} while (tries < 32 && field.power(weilPairing, subPrimePower).equals(field.one()));
					if (tries < 32) {
						result.add(torsionPoint2);
					}
				}
				torsionPointBasis.put(torsion, result);
			}
		}
		return torsionPointBasis.get(torsion);
	}

	public static class EllipticCurveExtensionResult<T extends Element<T>, B extends Element<B>, E extends AlgebraicExtensionElement<B, E>, FE extends FieldExtension<B, E, FE>> {
		private EllipticCurve<E> curve;
		private MathMap<ProjectivePoint<T>, ProjectivePoint<E>> embedding;

		private EllipticCurveExtensionResult(EllipticCurve<E> curve,
				MathMap<ProjectivePoint<T>, ProjectivePoint<E>> embedding) {
			this.curve = curve;
			this.embedding = embedding;
		}

		public EllipticCurve<E> getCurve() {
			return curve;
		}

		public MathMap<ProjectivePoint<T>, ProjectivePoint<E>> getEmbedding() {
			return embedding;
		}
	}

	public <B extends Element<B>, E extends AlgebraicExtensionElement<B, E>, FE extends FieldExtension<B, E, FE>> EllipticCurveExtensionResult<T, B, E, FE> extendBaseField(
			Extension<T, B, E, FE> extension) {
		MathMap<T, E> embedding = extension.embeddingMap();
		EllipticCurve<E> curve = new EllipticCurve<E>(extension.extension(), embedding.evaluate(a1),
				embedding.evaluate(a2), embedding.evaluate(a3), embedding.evaluate(a4), embedding.evaluate(a6));
		if (numberOfPoints != null) {
			curve.countPointsUpwards(numberOfPoints, field.getNumberOfElements());
		}
		return new EllipticCurveExtensionResult<>(curve, new MathMap<>() {
			@Override
			public ProjectivePoint<E> evaluate(ProjectivePoint<T> t) {
				return new ProjectivePoint<>(extension.extension(), embedding.evaluate(t.getCoord(1)),
						embedding.evaluate(t.getCoord(2)), embedding.evaluate(t.getCoord(3)));
			}
		});
	}

	public IntE embeddingDegree() {
		Integers z = Integers.z();
		FactorizationResult<IntE, IntE> factors = z.uniqueFactorization(z.getInteger(getNumberOfElements()));
		IntE embeddingDegree = z.one();
		for (IntE prime : factors.primeFactors()) {
			IntE torsion = z.power(prime, factors.multiplicity(prime));
			List<ProjectivePoint<T>> torsionPointBasis = getTorsionPointBasis(torsion);
			while (torsionPointBasis.isEmpty()) {
				torsion = z.divideChecked(torsion, prime);
				torsionPointBasis = getTorsionPointBasis(torsion);

			}
			embeddingDegree = z.lcm(embeddingDegree(torsion), embeddingDegree);
		}
		return embeddingDegree;
	}

	public IntE embeddingDegree(IntE torsion) {
		Integers z = Integers.z();
		if (z.isDivisible(torsion, z.getInteger(field.characteristic()))) {
			return z.one();
		}
		List<ProjectivePoint<T>> torisonPointBasis = getTorsionPointBasis(torsion);
		if (torisonPointBasis.size() == 0) {
			throw new ArithmeticException("No torsion point found!");
		}
		if (torisonPointBasis.size() == 2) {
			return z.one();
		}
		ModuloIntegerRing mod = new ModuloIntegerRing(torsion.getValue());
		IntE order = mod.getOrder(mod.getElement(field.getNumberOfElements()));
		if (z.gcd(torsion, z.getInteger(field.getNumberOfUnits())).equals(z.one())) {
			return order;
		}
		int multiplier = 0;
		ModuloIntegerRingElement sum = mod.zero();
		ModuloIntegerRingElement q = mod.getInteger(field.getNumberOfElements());
		do {
			for (int i = multiplier * order.intValueExact(); i < (multiplier + 1) * order.intValueExact(); i++) {
				sum = mod.add(mod.power(q, i), sum);
			}
			multiplier++;
		} while (!sum.equals(mod.zero()));
		return z.multiply(multiplier, order);
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
		if (minimalDefinedDegree < 0) {
			Integers z = Integers.z();
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
			for (IntE possibleDegree : z.factors(z.getInteger(degree))) {
				int i = possibleDegree.intValueExact();
				boolean degreeFound = true;
				for (T element : elements) {
					if (!field.power(element, p.pow(i)).equals(element)) {
						degreeFound = false;
						break;
					}
				}
				if (degreeFound) {
					minimalDefinedDegree = i;
					break;
				}
			}
		}
		return minimalDefinedDegree;
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
		if (domain.numberOfPoints != null) {
			this.numberOfPoints = domain.numberOfPoints;
		}
	}

	public void setNumberOfElements(BigInteger numberOfElements) {
		for (int i = 0; i < 10; i++) {
			if (!multiply(numberOfElements, getRandomElement()).equals(neutral())) {
				return;
			}
		}
		this.numberOfPoints = numberOfElements;
	}

	@Override
	public BigInteger getNumberOfElements() {
		if (!this.field.isFinite()) {
			return BigInteger.valueOf(-1);
		}
		if (this.numberOfPoints != null) {
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
			if ((l.intValue() - 1) / 2 < qbar) {
				qbar -= l.intValue();
			}
			List<Polynomial<T>> idealGen = new ArrayList<Polynomial<T>>();
			Polynomial<T> psiL = this.getDivisionPolynomial(l.intValue());
			idealGen.add(r.getEmbedding(psiL, new int[] { 0, 1 }));
			idealGen.add(r.getEmbedding(this.definingPolynomial, new int[] { 0, 1 }));
			PolynomialIdeal<T> ideal = r.getIdeal(idealGen);
			psiL = this.univariateRing.getEmbedding(psiL, new int[] { 0 });
			if (!field.roots(psiL).isEmpty()) {
				T root = field.roots(psiL).keySet().iterator().next();
				if (!getPossibleY(root).isEmpty()) {
					// P -t*P + qbar*P = 0
					// t = qbar + 1
					t.add(BigInteger.valueOf(qbar + 1));
					continue;
				}
			}
			CoordinateRing<T> cr = ideal.divideOut();
			// Calculating x^q, y^q, x^(q^2), y^(q^2)
			CoordinateRingElement<T> xq = cr.power(cr.getEmbedding(x), q);
			CoordinateRingElement<T> yq = cr.power(cr.getEmbedding(y), q);
			List<CoordinateRingElement<T>> xyq = new ArrayList<>();
			xyq.add(xq);
			xyq.add(yq);
			CoordinateRingElement<T> xqq = cr.power(xq, q);
			CoordinateRingElement<T> yqq = cr.power(yq, q);
			List<CoordinateRingElement<T>> qxy = multiplyGeneric(qbar, cr);
			CoordinateRingElement<T> qx = qxy.get(0);
			CoordinateRingElement<T> qy = qxy.get(1);
			if (!xqq.equals(qx)) {
				List<CoordinateRingElement<T>> lhs = this.addGeneric(xqq, yqq, qx, qy, cr);
				CoordinateRingElement<T> lhsX = lhs.get(0);
				CoordinateRingElement<T> lhsY = lhs.get(1);
				boolean found = false;
				for (int tc = 1; tc <= (l.intValue() - 1) / 2; tc++) {
					List<CoordinateRingElement<T>> candidate = this.multiplyGeneric(tc, cr);
					CoordinateRingElement<T> tcqx = cr.power(candidate.get(0), q);
					if (tcqx.equals(lhsX)) {
						found = true;
						CoordinateRingElement<T> tcqy = cr.power(candidate.get(1), q);
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
					List<CoordinateRingElement<T>> candidate = this.multiplyGeneric(w.getValue().intValueExact(), cr);
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

//	private Polynomial<T> invertXOnlyPcolynomial(Polynomial<T> a, Polynomial<T> psiL) {
//		PolynomialRing<T> univarring = this.univariateRing;
//		Polynomial<T> univarA = univarring.getEmbedding(a, new int[] { 0 });
//		// Inverting psiQsq mod psiL:
//		ExtendedEuclideanResult<Polynomial<T>> egcd = univarring.extendedEuclidean(univarA, psiL);
//		Polynomial<T> univarAInv = univarring.multiply(egcd.getCoeff1(), univarring.inverse(egcd.getGcd()));
//		return this.affineRing.getEmbedding(univarAInv, new int[] { 0 });
//	}

	private List<CoordinateRingElement<T>> addGeneric(CoordinateRingElement<T> x1, CoordinateRingElement<T> y1,
			CoordinateRingElement<T> x2, CoordinateRingElement<T> y2, CoordinateRing<T> cr) {
		CoordinateRingElement<T> ydiff = cr.subtract(y2, y1);
		CoordinateRingElement<T> xdiff = cr.subtract(x2, x1);
		CoordinateRingElement<T> s = cr.divideChecked(ydiff, xdiff);
		CoordinateRingElement<T> xr = cr.add(cr.multiply(s, s), cr.multiply(a1, s),
				cr.getEmbedding(field.multiply(-1, a2)));
		xr = cr.subtract(xr, cr.add(x1, x2));
		CoordinateRingElement<T> yr = cr.add(y1, cr.multiply(s, cr.subtract(xr, x1)));
		List<CoordinateRingElement<T>> result = new ArrayList<>();
		result.add(xr);
		result.add(cr.add(cr.negative(yr), cr.negative(cr.multiply(a1, xr)), cr.negative(cr.getEmbedding(a3))));
		return result;
	}

	private List<CoordinateRingElement<T>> multiplyGeneric(int n, CoordinateRing<T> cr) {
		// Calculating [qbar] (x, y)
		PolynomialRing<T> r = this.affineRing;
		int nabs = Math.abs(n);
		Polynomial<T> psiNsq = this.getDivisionPolynomial(nabs);
		psiNsq = r.multiply(psiNsq, psiNsq);
		psiNsq = cr.getPolynomialRing()
				.getEmbedding(r.reduce(psiNsq, Collections.singletonList(this.definingPolynomial)), new int[] { 0, 1 });
		CoordinateRingElement<T> psiNsqInv = cr.inverse(cr.getEmbedding(psiNsq));// .invertXOnlyPolynomial(psiNsq,
																					// psiL));
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

	private Polynomial<Ext> secondModularPolynomial(PolynomialRing<Ext> r) {
		Integers z = Integers.z();
		Polynomial<Ext> phi2 = r.add(r.getVarPower(1, 3), r.getVarPower(2, 3));
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

	private Polynomial<Ext> thirdModularPolynomial(PolynomialRing<Ext> r) {
		Polynomial<Ext> phi3 = r.add(r.getVarPower(1, 4), r.getVarPower(2, 4));
		// [1,0] 1855425871872000000000
		// [1,1] -770845966336000000
		// [2,0] 452984832000000
		// [2,1] 8900222976000
		// [2,2] 2587918086
		// [3,0] 36864000
		// [3,1] -1069956
		// [3,2] 2232
		// [3,3] -1
		// [4,0] 1
		phi3 = r.subtract(phi3, r.multiply(r.getVarPower(1, 3), r.getVarPower(2, 3)));
		phi3 = r.add(phi3, r.multiply(2232, r.add(r.multiply(r.getVarPower(1, 3), r.getVarPower(2, 2)),
				r.multiply(r.getVarPower(1, 2), r.getVarPower(2, 3)))));
		phi3 = r.subtract(phi3, r.multiply(1069956,
				r.add(r.multiply(r.getVarPower(1, 3), r.getVar(2)), r.multiply(r.getVar(1), r.getVarPower(2, 3)))));
		phi3 = r.add(phi3, r.multiply(36864000, r.add(r.getVarPower(1, 3), r.getVarPower(2, 3))));
		phi3 = r.add(phi3, r.multiply(BigInteger.valueOf(2587918086L), r.getVarPower(1, 2), r.getVarPower(2, 2)));
		phi3 = r.add(phi3, r.multiply(new BigInteger("8900222976000"),
				r.add(r.multiply(r.getVarPower(1, 2), r.getVar(2)), r.multiply(r.getVar(1), r.getVarPower(2, 2)))));
		phi3 = r.add(phi3,
				r.multiply(new BigInteger("452984832000000"), r.add(r.getVarPower(1, 2), r.getVarPower(2, 2))));
		phi3 = r.subtract(phi3, r.multiply(new BigInteger("770845966336000000"), r.getVar(1), r.getVar(2)));
		phi3 = r.add(phi3, r.multiply(new BigInteger("1855425871872000000000"), r.add(r.getVar(1), r.getVar(2))));
		return phi3;
	}

	// https://www.ams.org/journals/mcom/2003-72-241/S0025-5718-02-01434-5/S0025-5718-02-01434-5.pdf
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
		// We now have an elliptic curve with a3 = a4 = a6 = 0 and a1 = 1, i.e. only a2
		// differs between curves, and j(E) not in F4

		// Construction p-adic lift of the field:
		FiniteField f2q = (FiniteField) field;
		FFE j = (FFE) jInvariant();
		if (f2q.degree() != degree) {
			FiniteField jInvariantField = FiniteField.getFiniteField(f2q.minimalPolynomial(j),
					PrimeField.getPrimeField(2));
			countPointsUpwards(fromJInvariant(jInvariantField, jInvariantField.alpha()).getNumberOfElements(),
					jInvariantField.getNumberOfElements());
			return numberOfPoints;
		}
		FFE jPower = j;
		List<FFE> jInvariants = new ArrayList<>();
		do {
			jInvariants.add(jPower);
			jPower = f2q.power(jPower, 2);
		} while (!jPower.equals(j));
		int precision = MiscAlgorithms.DivRoundUp(degree + 3, 2) + 10;
		CompletedNumberField lift = getFieldLift2(precision);
		List<Ext> liftedJInvariants = liftJInvariants2(jInvariants, lift, precision);
		Ext uSqr = lift.one();
		for (int i = 0; i < liftedJInvariants.size(); i++) {
			Ext jI = liftedJInvariants.get(i);
			Ext jI1 = liftedJInvariants.get((i + 1) % liftedJInvariants.size());
			Ext xCoordinateNumerator = lift.add(lift.multiply(jI, jI), lift.multiply(195120, jI),
					lift.multiply(4095, jI1), lift.getInteger(660960000));
			xCoordinateNumerator = lift.multiply(-1, xCoordinateNumerator);
			Ext xCoordinateDenominator = lift.add(lift.multiply(jI, jI),
					lift.multiply(-1, jI1, lift.add(lift.multiply(512, jI), lift.getInteger(-372735))),
					lift.multiply(563760, jI), lift.getInteger(BigInteger.valueOf(8981280000L)));
			xCoordinateDenominator = lift.multiply(8, xCoordinateDenominator);
			Ext xCoordinateHalf = lift.divide(xCoordinateNumerator, xCoordinateDenominator);
			Ext yCoordinate = lift.negative(xCoordinateHalf);
			Ext t = lift.subtract(lift.multiply(12, xCoordinateHalf, xCoordinateHalf),
					lift.add(lift.divide(lift.getInteger(36), lift.subtract(jI1, lift.getInteger(1728))), yCoordinate));
			Ext a = lift.add(lift.divide(lift.getInteger(-36), lift.subtract(jI1, lift.getInteger(1728))),
					lift.multiply(-5, t));
			Ext b = lift.add(lift.negative(lift.inverse(lift.subtract(jI1, lift.getInteger(1728)))),
					lift.multiply(-1, lift.add(lift.one(), lift.multiply(14, xCoordinateHalf)), t));
			Ext uSqrI = lift.negative(lift.divide(lift.subtract(lift.multiply(48, a), lift.one()),
					lift.add(lift.multiply(6 * 12 * 12, b), lift.multiply(-6 * 12, a), lift.one())));
			uSqr = lift.multiply(uSqrI, uSqr);
		}
		Ext cSqr = lift.round(lift.inverse(uSqr), MiscAlgorithms.DivRoundUp(degree + 3, 2) + 1);
		Map<PAdicNumber, Integer> candidates = lift.getBaseField().sqrt(lift.asBaseFieldElement(cSqr));
		if (candidates.size() != 2) {
			throw new ArithmeticException("Could not find square root!");
		}
		Integers z = Integers.z();
		for (PAdicNumber c : candidates.keySet()) {
			IntE trace = lift.getBaseField().roundToInteger(c, MiscAlgorithms.DivRoundUp(degree + 3, 2));
			if (z.remainder(trace, z.getInteger(4)).equals(z.one())) {
				numberOfPoints = z.subtract(z.add(z.getInteger(field.getNumberOfElements()), z.one()), trace)
						.getValue();
				return numberOfPoints;
			}
		}
		throw new ArithmeticException("Could not find candidate!");
	}

	private List<Ext> liftJInvariants2(List<FFE> reduced, CompletedNumberField liftedField, int desiredPrecision) {
		PolynomialRing<Ext> r = AbstractPolynomialRing.getPolynomialRing(liftedField, reduced.size(), Monomial.GREVLEX);
		List<Polynomial<Ext>> modularFunction = new ArrayList<>();
		Polynomial<Ext> modularPolynomial = secondModularPolynomial(
				AbstractPolynomialRing.getPolynomialRing(liftedField, 2, Monomial.GREVLEX));
		for (int i = 0; i < reduced.size(); i++) {
			List<Polynomial<Ext>> substitute = new ArrayList<>();
			substitute.add(r.getVar(i + 1));
			substitute.add(r.getVar((i + 1) % reduced.size() + 1));
			modularFunction.add(r.substitute(modularPolynomial, substitute));
		}
		return liftedField.ringOfIntegers()
				.henselLiftVector(r, new Vector<>(modularFunction), new Vector<>(reduced), desiredPrecision).asList();
	}

	private CompletedNumberField getFieldLift2(int accuracy) {
		FiniteField f2q = (FiniteField) field;
		UnivariatePolynomial<PFE> minimalPolynomial = f2q.minimalPolynomial();
		PAdicField q2 = new PAdicField(2, accuracy);
		UnivariatePolynomial<PAdicNumber> liftedMinimalPolynomial = q2.ringOfIntegers()
				.liftUnivariatePolynomial(minimalPolynomial);
		return q2.getExtension(liftedMinimalPolynomial).extension();
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
		if (!a2.equals(field.one()) && field.hasSqrt(a2)) {
			T sqrt = field.sqrt(a2).keySet().iterator().next();
			Isomorphism<T> isomorphism = new Isomorphism<>(this, sqrt, field.zero(), field.zero(), field.zero());
			numberOfPoints = isomorphism.getRange().getNumberOfElements();
			return numberOfPoints;
		} else if (!a2.equals(field.one())) {
			numberOfPoints = field.getNumberOfElements().add(BigInteger.ONE).multiply(BigInteger.TWO)
					.subtract(getQuadraticTwist().getNumberOfElements());
			return numberOfPoints;
		}
		if (degree == 1) {
			int count = 1;
			for (int i = 0; i < 3; i++) {
				T x = field.getInteger(i);
				List<T> ys = getPossibleY(x);
				for (T y : ys) {
					if (field.power(y, 3).equals(y)) {
						count++;
					}
				}
			}
			countPointsUpwards(BigInteger.valueOf(count), BigInteger.valueOf(3));
			return numberOfPoints;
		}
		// We now have an elliptic curve with a1 = a3 = a4 = 0 and a2 = 1, i.e. only a6
		// differs between curves, and j(E) not in F3

		// Construction p-adic lift of the field:
		FiniteField f3q = (FiniteField) field;
		FFE j = (FFE) jInvariant();
		if (f3q.degree() != degree) {
			FiniteField jInvariantField = FiniteField.getFiniteField(f3q.minimalPolynomial(j),
					PrimeField.getPrimeField(3));
			countPointsUpwards(fromJInvariant(jInvariantField, jInvariantField.alpha()).getNumberOfElements(),
					jInvariantField.getNumberOfElements());
			return numberOfPoints;
		}
		FFE jPower = j;
		List<FFE> jInvariants = new ArrayList<>();
		do {
			jInvariants.add(jPower);
			jPower = f3q.power(jPower, 3);
		} while (!jPower.equals(j));
		int precision = MiscAlgorithms.DivRoundUp(degree + 3, 2) + 10;
		CompletedNumberField lift = getFieldLift3(precision);
		List<Ext> liftedJInvariants = liftJInvariants3(jInvariants, lift, precision);
		Ext prevCoeff = lift
				.inverse(lift.subtract(liftedJInvariants.get(liftedJInvariants.size() - 1), lift.getInteger(1728)));
		EllipticCurve<Ext> prevCurve = new EllipticCurve<>(lift, lift.zero(), lift.one(), lift.zero(),
				lift.multiply(-576, prevCoeff), lift.multiply(-64, prevCoeff));
		Ext uSqr = lift.one();
		for (Ext liftedJ : liftedJInvariants) {
			Ext coeff = lift.inverse(lift.subtract(liftedJ, lift.getInteger(1728)));
			EllipticCurve<Ext> liftedCurve = new EllipticCurve<>(lift, lift.zero(), lift.one(), lift.zero(),
					lift.multiply(-576, coeff), lift.multiply(-64, coeff));
			for (ProjectivePoint<Ext> torsionPoint : liftedCurve.getTorsionPoints(3)) {
				if (torsionPoint.equals(liftedCurve.neutral())) {
					continue;
				}
				Value v = Value.INFINITY;
				for (Ext p : torsionPoint.getCoords()) {
					v = v.min(lift.valuation(p));
				}
				if (v.compareTo(Value.ZERO) < 0) {
					continue;
				}
				Ext pointX = torsionPoint.getDehomogenisedCoord(1, 3);
				Ext pointY = torsionPoint.getDehomogenisedCoord(2, 3);
				Ext gxp = lift.add(lift.multiply(3, pointX, pointX), lift.multiply(2, pointX), liftedCurve.getA4());
				Ext gyp = lift.multiply(-2, pointY);
				Ext t = lift.multiply(2, gxp);
				Ext up = lift.multiply(gyp, gyp);
				Ext u = lift.add(up, lift.multiply(pointX, t));
				Ext rangeA4 = lift.subtract(liftedCurve.getA4(), lift.multiply(5, t));
				Ext rangeA6 = lift.subtract(liftedCurve.getA6(),
						lift.add(lift.multiply(7, u), lift.multiply(liftedCurve.getB2(), t)));
				Ext targetA4 = prevCurve.getA4();
				Ext targetA6 = prevCurve.getA6();
				// x0 = (27*a4*a6Target + -1*a4 + -27*a6*a4Target + 9*a6 + a4Target +
				// -9*a6Target)/(27*a4*a4Target + -81*a4*a6Target + -6*a4 + -9*a4Target +
				// 27*a6Target + 2)
				Ext xOffsetNumerator = lift.add(
						lift.add(lift.multiply(27, rangeA4, targetA6), lift.multiply(-27, targetA4, rangeA4),
								lift.multiply(-1, rangeA4)),
						lift.add(targetA4, lift.multiply(9, rangeA6), lift.multiply(-9, targetA6)));
				Ext xOffsetDenominator = lift.add(
						lift.add(lift.multiply(27, rangeA4, targetA4), lift.multiply(-81, rangeA4, targetA6),
								lift.multiply(-6, rangeA4)),
						lift.add(lift.multiply(-9, targetA4), lift.multiply(27, targetA6), lift.getInteger(2)));
				Ext xOffset = lift.divide(xOffsetNumerator, xOffsetDenominator);
				uSqr = lift.multiply(lift.add(lift.multiply(3, xOffset), lift.one()), uSqr);
				break;
			}
			prevCurve = liftedCurve;
		}
		PAdicNumber cSqr = lift.asBaseFieldElement(lift.round(uSqr, MiscAlgorithms.DivRoundUp(degree + 3, 2) + 1));
		for (PAdicNumber c : lift.getBaseField().sqrt(cSqr).keySet()) {
			IntE trace = lift.getBaseField().roundToInteger(c, MiscAlgorithms.DivRoundUp(degree + 3, 2));
			Integers z = Integers.z();
			if (!z.remainder(trace, z.getInteger(3)).equals(z.one())) {
				continue;
			}
			numberOfPoints = z.subtract(z.add(z.getInteger(field.getNumberOfElements()), z.one()), trace).getValue();
			return numberOfPoints;
		}
		throw new ArithmeticException("No square root found!");
	}

	private List<Ext> liftJInvariants3(List<FFE> reduced, CompletedNumberField liftedField, int desiredPrecision) {
		PolynomialRing<Ext> r = AbstractPolynomialRing.getPolynomialRing(liftedField, reduced.size(), Monomial.GREVLEX);
		List<Polynomial<Ext>> modularFunction = new ArrayList<>();
		Polynomial<Ext> modularPolynomial = thirdModularPolynomial(
				AbstractPolynomialRing.getPolynomialRing(liftedField, 2, Monomial.GREVLEX));
		for (int i = 0; i < reduced.size(); i++) {
			List<Polynomial<Ext>> substitute = new ArrayList<>();
			substitute.add(r.getVar(i + 1));
			substitute.add(r.getVar((i + 1) % reduced.size() + 1));
			modularFunction.add(r.substitute(modularPolynomial, substitute));
		}
		return liftedField.ringOfIntegers()
				.henselLiftVector(r, new Vector<>(modularFunction), new Vector<>(reduced), desiredPrecision).asList();
	}

	private CompletedNumberField getFieldLift3(int accuracy) {
		FiniteField f3q = (FiniteField) field;
		UnivariatePolynomial<PFE> minimalPolynomial = f3q.minimalPolynomial();
		PAdicField q3 = new PAdicField(3, accuracy);
		UnivariatePolynomial<PAdicNumber> liftedMinimalPolynomial = q3.ringOfIntegers()
				.liftUnivariatePolynomial(minimalPolynomial);
		return q3.getExtension(liftedMinimalPolynomial).extension();
	}

	public ProjectiveMorphism<T> translationMorphism(ProjectivePoint<T> point) {
		if (point.equals(pointAtInfinity)) {
			List<Polynomial<T>> values = new ArrayList<>();
			values.add(projectiveRing.getVar(1));
			values.add(projectiveRing.getVar(2));
			values.add(projectiveRing.getVar(3));
			return new ProjectiveMorphism<>(asGenericProjectiveScheme(), asGenericProjectiveScheme(), values);
		}
		FunctionField<T> ff = getFunctionField();
		CoordinateRing<T> cr = getCoordinateRing();
		T x = point.getDehomogenisedCoord(1, 3);
		T y = point.getDehomogenisedCoord(2, 3);
		RationalFunction<T> s = ff.getFunction(cr.subtract(cr.getVar(2), cr.getEmbedding(y)),
				cr.subtract(cr.getVar(1), cr.getEmbedding(x)));
		RationalFunction<T> xResult = ff.add(ff.add(ff.multiply(s, s), ff.negative(ff.getEmbedding(a2))),
				ff.add(ff.multiply(ff.getEmbedding(a1), s), ff.negative(ff.getEmbedding(x)),
						ff.negative(ff.getFunction(cr.getVar(1), cr.one()))));
		RationalFunction<T> yThird = ff.add(ff.multiply(s, ff.subtract(xResult, ff.getEmbedding(x))),
				ff.getEmbedding(y));
		RationalFunction<T> yResult = ff
				.negative(ff.add(yThird, ff.multiply(ff.getEmbedding(a1), xResult), ff.getEmbedding(a3)));
		List<RationalFunction<T>> asRationalFunctions = new ArrayList<>();
		asRationalFunctions.add(xResult);
		asRationalFunctions.add(yResult);
		return ProjectiveMorphism.fromRationalFunctions(this.asGenericProjectiveScheme(),
				this.asGenericProjectiveScheme(), asRationalFunctions, 2);
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
		if (!a1.equals(field.zero()) || !a3.equals(field.zero())) {
			T minusHalf = field.inverse(field.getInteger(-2));
			T xyOffset = field.multiply(minusHalf, a1);
			T yOffset = field.multiply(minusHalf, a3);
			return new Isomorphism<>(this, field.one(), field.zero(), xyOffset, yOffset).getRange().getAdjustedRhs();
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
				for (T y : getPossibleY(x)) {
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
		return "Y^2" + coefficientToString(a1, "XY") + coefficientToString(a3, "Y") + " = X^3"
				+ coefficientToString(a2, "X^2") + coefficientToString(a4, "X") + coefficientToString(a6, "");
	}

	private String coefficientToString(T coeff, String monomial) {
		if (coeff.equals(field.zero())) {
			return "";
		}
		String toString = coeff.toString();
		if (toString.contains(" ")) {
			toString = "(" + toString + ")";
		}
		return " + " + toString + (monomial.length() > 0 ? "*" + monomial : "");
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
	public List<ProjectiveMorphism<T>> irreducibleComponents() {
		return Collections.singletonList(identityMorphism());
	}

	@Override
	public ProjectiveMorphism<T> reduced() {
		return identityMorphism();
	}
}
