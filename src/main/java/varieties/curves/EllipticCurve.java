package varieties.curves;

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import fields.finitefields.FiniteField;
import fields.finitefields.FiniteField.FFE;
import fields.finitefields.PrimeField;
import fields.finitefields.PrimeField.PFE;
import fields.floatingpoint.Complex;
import fields.floatingpoint.Complex.ComplexNumber;
import fields.helper.CoordinateRing;
import fields.helper.CoordinateRing.CoordinateRingElement;
import fields.helper.FieldEmbedding;
import fields.integers.Integers;
import fields.integers.Integers.IntE;
import fields.interfaces.Element;
import fields.interfaces.Field;
import fields.interfaces.Group;
import fields.interfaces.MathMap;
import fields.interfaces.Polynomial;
import fields.interfaces.PolynomialRing;
import fields.interfaces.Ring.ExtendedEuclideanResult;
import fields.interfaces.Ring.QuotientAndRemainderResult;
import fields.interfaces.UnivariatePolynomial;
import fields.interfaces.UnivariatePolynomialRing;
import fields.polynomials.AbstractPolynomialRing;
import fields.polynomials.Monomial;
import fields.polynomials.PolynomialIdeal;
import fields.vectors.Vector;
import util.ConCatMap;
import util.MiscAlgorithms;
import varieties.FunctionField;
import varieties.ProjectiveMorphism;
import varieties.ProjectivePoint;
import varieties.ProjectiveVariety;
import varieties.RationalFunction;
import varieties.curves.DivisorGroup.Divisor;
import varieties.projective.AbstractProjectiveVariety;

public class EllipticCurve<T extends Element<T>> extends AbstractProjectiveVariety<T>
		implements SmoothCurve<T>, Group<ProjectivePoint<T>> {
	private static class Isomorphism<T extends Element<T>> implements Isogeny<T> {
		private final EllipticCurve<T> domain;
		private final EllipticCurve<T> range;
		private final Field<T> field;
		private final T u;

		private Isomorphism(EllipticCurve<T> domain, T u) {
			this.u = u;
			this.domain = domain;
			this.field = domain.getField();
			this.range = new EllipticCurve<T>(field, field.multiply(field.power(u, 4), domain.getA()),
					field.multiply(field.power(u, 6), domain.getB()));
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
			T x = point.getCoord(1);
			T y = point.getCoord(2);
			T z = point.getCoord(3);
			return new ProjectivePoint<>(field, field.multiply(field.power(u, 2), x),
					field.multiply(field.power(u, 3), y), z);
		}

		@Override
		public ProjectiveMorphism<T> asMorphism() {
			PolynomialRing<T> projectiveRing = domain.asProjectiveVariety().homogenousPolynomialRing();
			List<Polynomial<T>> asPolynomials = new ArrayList<>();
			asPolynomials.add(projectiveRing.multiply(field.power(u, 2), projectiveRing.getVar(1)));
			asPolynomials.add(projectiveRing.multiply(field.power(u, 3), projectiveRing.getVar(2)));
			asPolynomials.add(projectiveRing.getVar(3));
			return new ProjectiveMorphism<>(domain.asProjectiveVariety(), range.asProjectiveVariety(), asPolynomials);
		}

		@Override
		public BigInteger getDegree() {
			return BigInteger.ONE;
		}

		@Override
		public Isogeny<T> getDual() {
			return new Isomorphism<>(range, field.inverse(u));
		}
	}

	private Field<T> field;
	private T a;
	private T b;
	private ProjectivePoint<T> pointAtInfinity;
	private PolynomialRing<T> projectiveRing;
	private PolynomialRing<T> affineRing;
	private UnivariatePolynomialRing<T> univariateRing;
	private CoordinateRing<T> coordinateRing;
	private FunctionField<T> functionField;
	private Polynomial<T> rhs;
	private Map<Integer, Polynomial<T>> divisionPolynomials;
	private Polynomial<T> definingPolynomial;
	private BigInteger numberOfPoints = BigInteger.valueOf(-1);
	private boolean superSingular;
	private boolean superSingularDeterminied;

	private static <T extends Element<T>> ProjectiveVariety<T> asProjectiveVariety(Field<T> field, T a, T b) {
		PolynomialRing<T> ring = AbstractPolynomialRing.getPolynomialRing(field, 3, Monomial.GREVLEX);
		Polynomial<T> f = ring.getEmbedding(field.one(), new int[] { 3, 0, 0 });
		f = ring.add(f, ring.getEmbedding(field.negative(field.one()), new int[] { 0, 2, 1 }));
		f = ring.add(f, ring.getEmbedding(a, new int[] { 1, 0, 2 }));
		f = ring.add(f, ring.getEmbedding(b, new int[] { 0, 0, 3 }));
		return new ProjectiveVariety<>(field, ring, Collections.singletonList(f));
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

	public static FromPolynomialResult<FFE> fromPolynomial(FiniteField field, Polynomial<FFE> polynomial) {
		PolynomialRing<FFE> ring = AbstractPolynomialRing.getPolynomialRing(field, 2, Monomial.REVLEX);
		if (polynomial.numberOfVariables() == 3 && polynomial.getPolynomialRing().isHomogeneous(polynomial)) {
			polynomial = ring.getEmbedding(polynomial.getPolynomialRing().dehomogenize(polynomial, 3),
					new int[] { 0, 1 });
		} else if (polynomial.numberOfVariables() == 2) {
			polynomial = ring.getEmbedding(polynomial, new int[] { 0, 1 });
		} else {
			throw new RuntimeException("Unexpected number of variables");
		}
		FieldEmbedding<PFE, FFE, FiniteField> embedding = new FieldEmbedding<>(field);
		BigInteger ch = field.characteristic();
		Polynomial<FFE> pPower = polynomial;
		int power = 0;
		do {
			power++;
			Map<Monomial, FFE> coeffs = new TreeMap<>();
			for (Monomial m : pPower.monomials()) {
				coeffs.put(m, field.power(pPower.coefficient(m), ch));
			}
			pPower = pPower.getPolynomialRing().getPolynomial(coeffs);
		} while (!pPower.equals(polynomial));
		if (!field.getNumberOfElements().equals(ch.pow(power))) {
			FiniteField embed = FiniteField.getFiniteField(ch, power);
			embedding = new FieldEmbedding<>(field, embed);
			ring = AbstractPolynomialRing.getPolynomialRing(embed, 2, Monomial.REVLEX);
			polynomial = ring.getEmbedding(polynomial,
					new ConCatMap<>(embedding.asVectorMap(), new MathMap<Vector<FFE>, FFE>() {
						@Override
						public FFE evaluate(Vector<FFE> t) {
							return t.get(1);
						}
					}));
		}
		FiniteField f = embedding.getEmbeddedField();
		FFE y2 = f.zero();
		FFE x3 = f.zero();
		for (Monomial m : polynomial.monomials()) {
			if (polynomial.coefficient(m).equals(f.zero())) {
				continue;
			}
			if (2 * m.exponents()[0] + 3 * m.exponents()[1] > 6) {
				throw new RuntimeException("degrees wrong");
			}
			if (m.exponents()[0] == 3) {
				x3 = polynomial.coefficient(m);
			}
			if (m.exponents()[1] == 2) {
				y2 = polynomial.coefficient(m);
			}
		}
		if (x3.equals(f.zero()) || y2.equals(f.zero())) {
			throw new RuntimeException("No X^3 or Y^2 term!");
		}
		y2 = f.divide(y2, x3);
		while (!f.hasSqrt(y2)) {
			FFE rng = f.power(f.getRandomElement(), 3);
			if (rng.equals(f.zero())) {
				continue;
			}
			y2 = f.multiply(y2, rng);
		}
		if (!y2.equals(polynomial.leadingCoefficient())) {
			polynomial = ring.scalarMultiply(f.divide(y2, polynomial.leadingCoefficient()), polynomial);
		}
		y2 = polynomial.coefficient(ring.getMonomial(new int[] { 0, 2 }));
		x3 = polynomial.coefficient(ring.getMonomial(new int[] { 3, 0 }));
		List<Polynomial<FFE>> substitutions = new ArrayList<>();
		List<Polynomial<FFE>> sub = new ArrayList<>();
		sub.add(ring.getEmbedding(f.negative(f.inverse(f.roots(x3, 3).keySet().iterator().next())),
				new int[] { 1, 0 }));
		sub.add(ring.getEmbedding(f.inverse(f.roots(y2, 2).keySet().iterator().next()), new int[] { 0, 1 }));
		substitutions.addAll(sub);
		polynomial = ring.substitute(polynomial, sub);
		FFE xy = polynomial.coefficient(ring.getMonomial(new int[] { 1, 1 }));
		FFE y = polynomial.coefficient(ring.getMonomial(new int[] { 0, 1 }));
		xy = f.divide(xy, f.getInteger(2));
		y = f.divide(y, f.getInteger(2));
		Polynomial<FFE> ysub = ring.getVar(2);
		ysub = ring.subtract(ysub, ring.getEmbedding(xy, new int[] { 1, 0 }));
		ysub = ring.subtract(ysub, ring.getEmbedding(y, new int[] { 0, 0 }));
		sub.clear();
		sub.add(ring.getVar(1));
		sub.add(ysub);
		substitutions.set(1, ring.substitute(substitutions.get(1), sub));
		polynomial = ring.substitute(polynomial, sub);
		Polynomial<FFE> xsub = ring.getVar(1);
		xsub = ring.add(xsub, ring
				.getEmbedding(f.divide(polynomial.coefficient(ring.getMonomial(new int[] { 2, 0 })), f.getInteger(3))));
		sub.clear();
		sub.add(xsub);
		sub.add(ring.getVar(2));
		substitutions.set(0, ring.substitute(substitutions.get(0), sub));
		polynomial = ring.substitute(polynomial, sub);
		return new FromPolynomialResult<>(new EllipticCurve<>(field,
				field.negative(embedding.getEmbedding(polynomial.coefficient(ring.getMonomial(new int[] { 1, 0 })))),
				field.negative(embedding.getEmbedding(polynomial.coefficient(ring.getMonomial(new int[] { 0, 0 }))))),
				substitutions);
	}

	public static <T extends Element<T>> EllipticCurve<T> fromJInvariant(Field<T> field, T j) {
		if (j.equals(field.zero())) {
			return new EllipticCurve<>(field, field.zero(), field.one());
		}
		if (j.equals(field.getInteger(1728))) {
			return new EllipticCurve<>(field, field.negative(field.one()), field.zero());
		}
		T a = field.divide(field.multiply(4, field.subtract(field.getInteger(1728), j)), field.multiply(27, j));
		T b = field.power(a, 2);
		return new EllipticCurve<>(field, a, b);
	}

	public static <T extends Element<T>> EllipticCurve<T> fromLengendre(Field<T> field, T lambda) {
		T numerator = field.multiply(256,
				field.power(field.add(field.multiply(lambda, lambda), field.multiply(-1, lambda), field.one()), 3));
		T denominator = field.multiply(field.power(lambda, 2), field.power(field.subtract(lambda, field.one()), 2));
		T j = field.divide(numerator, denominator);
		return EllipticCurve.fromJInvariant(field, j);
//		UnivariatePolynomialRing<T> r = field.getUnivariatePolynomialRing();
//		UnivariatePolynomial<T> rhs = r.getPolynomial(field.zero(), field.negative(field.one()), field.one());
//		rhs = r.multiply(rhs, r.getPolynomial(field.negative(lambda), field.one()));
//		UnivariatePolynomial<T> sub = r.getPolynomial(
//				field.negative(field.divide(rhs.univariateCoefficient(2), field.getInteger(3))), field.one());
//		UnivariatePolynomial<T> newRhs = r.substitute(rhs, Collections.singletonList(sub));
//		return new EllipticCurve<T>(field, newRhs.univariateCoefficient(1), newRhs.univariateCoefficient(0));
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

	private static ComplexNumber capitalG(Complex c, ComplexNumber tau, int weight) {
		ComplexNumber result = c.zero();
		ComplexNumber prev;
		int i = 1;
		do {
			prev = result;
			for (int n = 0; n <= i; n++) {
				int m = i - n;
				ComplexNumber summand1 = c.power(c.inverse(c.add(c.getInteger(m), c.multiply(n, tau))), 2 * weight);
				ComplexNumber summand2 = c.power(c.inverse(c.add(c.getInteger(m), c.multiply(-n, tau))), 2 * weight);
				result = c.add(result, summand1, summand2);
				result = c.add(result, summand1, summand2);
			}
		} while (!c.close(prev, result));
		return result;
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
		if (tau.complexPart().compareTo(c.getBaseField().zero()) <= 0) {
			throw new ArithmeticException("Tau not in upper half plane!");
		}
		ComplexNumber g4 = c.multiply(60, capitalG(c, tau, 2));
		ComplexNumber g6 = c.multiply(140, capitalG(c, tau, 3));
		ComplexNumber j = c.multiply(1728,
				c.divide(c.power(g4, 3), c.subtract(c.power(g4, 3), c.multiply(27, g6, g6))));
		return new FromTauResult(tau, fromJInvariant(c, j), weierstrassP(c, tau), weierstrassPDerivative(c, tau));
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
		super(asProjectiveVariety(field, a, b));
		this.field = field;
		this.a = a;
		this.b = b;
		this.pointAtInfinity = new ProjectivePoint<T>(this.field, this.field.zero(), this.field.one(),
				this.field.zero());
		if (this.field
				.add(this.field.multiply(4, this.field.power(a, 3)), this.field.multiply(27, this.field.power(b, 2)))
				.equals(this.field.zero()))
			throw new ArithmeticException("Singular curve");
		this.projectiveRing = asProjectiveVariety().homogenousPolynomialRing();
		this.affineRing = AbstractPolynomialRing.getPolynomialRing(this.field, 2, Monomial.REVLEX);
		PolynomialRing<T> r = this.affineRing;
		Polynomial<T> x = r.getVar(1);
		Polynomial<T> x3 = r.getVarPower(1, 3);
		Polynomial<T> y2 = r.getVarPower(2, 2);
		Polynomial<T> rb = r.getEmbedding(this.b);
		Polynomial<T> rhs = r.add(x3, r.multiply(this.a, x), rb);
		this.functionField = new FunctionField<>(this);
		this.definingPolynomial = r.subtract(rhs, y2);
		this.univariateRing = this.field.getUnivariatePolynomialRing();
		r = this.univariateRing;
		x = r.getVar(1);
		x3 = r.getVarPower(1, 3);
		rb = r.getEmbedding(this.b);
		this.rhs = r.add(x3, r.multiply(this.a, x), rb);
	}

	public T jInvariant() {
		T fourACube = this.field.multiply(4, this.field.power(a, 3));
		T twentySevenBSquare = this.field.multiply(27, this.field.power(b, 2));
		return this.field.multiply(1728, this.field.divide(fourACube, this.field.add(fourACube, twentySevenBSquare)));
	}

	public List<Isogeny<T>> getIsomorphisms(EllipticCurve<T> range) {
		if (!this.jInvariant().equals(range.jInvariant())) {
			throw new ArithmeticException("Not isomorphic");
		}
		UnivariatePolynomial<T> toSolve;
		if (a.equals(field.zero())) {
			toSolve = univariateRing.toUnivariate(univariateRing.subtract(univariateRing.getVarPower(6),
					univariateRing.getEmbedding(field.divide(range.b, b))));
		} else if (b.equals(field.zero())) {
			toSolve = univariateRing.toUnivariate(univariateRing.subtract(univariateRing.getVarPower(4),
					univariateRing.getEmbedding(field.divide(range.a, a))));
		} else {
			toSolve = univariateRing.toUnivariate(univariateRing.subtract(univariateRing.getVarPower(2),
					univariateRing.getEmbedding(field.divide(field.multiply(range.b, a), field.multiply(b, range.a)))));
		}
		List<Isogeny<T>> isomorphisms = new ArrayList<>();
		for (T root : field.roots(toSolve).keySet()) {
			isomorphisms.add(new Isomorphism<>(this, root));
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
		if (q.equals(BigInteger.ONE)) {
			return frobenius;
		}
		List<Isogeny<T>> isomorphisms = frobenius.getRange().getIsomorphisms(this);
		return new CompositionIsogeny<>(frobenius, isomorphisms.get(0));
	}

	public BigInteger trace() {
		return field.getNumberOfElements().add(BigInteger.ONE).subtract(getNumberOfElements());
	}

	public boolean isSupersingular() {
		if (!superSingularDeterminied) {
			this.superSingular = attemptSupersingularCount();
			this.superSingularDeterminied = true;
		}
		return superSingular;
	}

	public T getA() {
		return a;
	}

	public T getB() {
		return b;
	}

	public CoordinateRing<T> getCoordinateRing() {
		if (this.coordinateRing == null) {
			this.coordinateRing = new CoordinateRing<T>(affineRing,
					affineRing.getIdeal(Collections.singletonList(definingPolynomial)));
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
		if (this.divisionPolynomials.containsKey(m)) {
			Polynomial<T> p = this.divisionPolynomials.get(m);
			if (m % 2 == 0) {
				return this.affineRing.multiply(2, this.affineRing.getVar(2), p);
			} else {
				return p;
			}
		}
		PolynomialRing<T> r = this.affineRing;
		Polynomial<T> o = r.one();
		Polynomial<T> x = r.getVar(1);
		Polynomial<T> x2 = r.getVarPower(1, 2);
		Polynomial<T> x3 = r.getVarPower(1, 3);
		Polynomial<T> x4 = r.getVarPower(1, 4);
		Polynomial<T> x6 = r.getVarPower(1, 6);
		Polynomial<T> a = r.getEmbedding(this.a);
		Polynomial<T> a2 = r.multiply(a, a);
		Polynomial<T> a3 = r.multiply(a2, a);
		Polynomial<T> b = r.getEmbedding(this.b);
		Polynomial<T> b2 = r.multiply(b, b);
		Polynomial<T> ab = r.multiply(a, b);
		Polynomial<T> y2 = r.add(x3, r.multiply(a, x), b);
		if (m == 0) {
			this.divisionPolynomials.put(0, r.zero());
		} else if (m == 1) {
			this.divisionPolynomials.put(1, o);
		} else if (m == 2) {
			this.divisionPolynomials.put(2, o);
		} else if (m == 3) {
			Polynomial<T> p = r.multiply(3, x4);
			p = r.add(p, r.multiply(6, a, x2));
			p = r.add(p, r.multiply(12, b, x));
			p = r.add(p, r.multiply(-1, a2));
			this.divisionPolynomials.put(3, p);
		} else if (m == 4) {
			Polynomial<T> p = x6;
			p = r.add(p, r.multiply(5, a, x4));
			p = r.add(p, r.multiply(20, b, x3));
			p = r.add(p, r.multiply(-5, a2, x2));
			p = r.add(p, r.multiply(-4, ab, x));
			p = r.add(p, r.multiply(-8, b2));
			p = r.add(p, r.multiply(-1, a3));
			p = r.multiply(2, p);
			this.divisionPolynomials.put(4, p);
		} else if (m % 2 == 0) {
			int n = m / 2;
			for (int i = n - 2; i <= n + 2; i++) {
				this.getDivisionPolynomial(i);
			}
			Polynomial<T> psiNm2 = this.divisionPolynomials.get(n - 2);
			Polynomial<T> psiNm1 = this.divisionPolynomials.get(n - 1);
			Polynomial<T> psiN = this.divisionPolynomials.get(n);
			Polynomial<T> psiN1 = this.divisionPolynomials.get(n + 1);
			Polynomial<T> psiN2 = this.divisionPolynomials.get(n + 2);
			Polynomial<T> p = r.subtract(r.multiply(psiN2, r.power(psiNm1, 2)), r.multiply(psiNm2, r.power(psiN1, 2)));
			p = r.multiply(psiN, p);
			this.divisionPolynomials.put(m, p);
		} else if (m % 2 == 1) {
			int n = (m - 1) / 2;
			for (int i = n - 1; i <= n + 2; i++) {
				this.getDivisionPolynomial(i);
			}
			Polynomial<T> psiNm1 = this.divisionPolynomials.get(n - 1);
			Polynomial<T> psiN = this.divisionPolynomials.get(n);
			Polynomial<T> psiN1 = this.divisionPolynomials.get(n + 1);
			Polynomial<T> psiN2 = this.divisionPolynomials.get(n + 2);
			Polynomial<T> firstTerm = r.multiply(psiN2, r.power(psiN, 3));
			Polynomial<T> secondTerm = r.multiply(psiNm1, r.power(psiN1, 3));
			if (n % 2 == 0) {
				firstTerm = r.multiply(16, firstTerm, y2, y2);
			} else {
				secondTerm = r.multiply(16, secondTerm, y2, y2);
			}
			this.divisionPolynomials.put(m, r.subtract(firstTerm, secondTerm));
		}
		return this.getDivisionPolynomial(m);
	}

	public ProjectivePoint<T> neutral() {
		return this.pointAtInfinity;
	}

	public ProjectivePoint<T> negative(ProjectivePoint<T> p) {
		return new ProjectivePoint<T>(this.field, p.getCoord(1), this.field.negative(p.getCoord(2)), p.getCoord(3));
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
			s = field.divide(field.add(field.multiply(3, px, px), a), field.multiply(2, py));
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
		return a.equals(ec.a) && b.equals(ec.b);
	}

	@Override
	public int getEmbeddingDimension() {
		return 2;
	}

	@Override
	public int dimension() {
		return 1;
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

	public RationalFunction<T> getRationalFunction(ProjectivePoint<T> zero,
			ProjectivePoint<T> pole1, ProjectivePoint<T> pole2) {
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

	public FunctionField<T> getFunctionField() {
		return functionField;
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
		if (px.equals(qx))
			s = this.field.divide(this.field.add(this.field.multiply(3, px, px), a), this.field.multiply(2, py));
		else
			s = this.field.divide(this.field.subtract(py, qy), this.field.subtract(px, qx));
		r1 = this.field.multiply(s, s);
		r1 = this.field.subtract(this.field.subtract(r1, px), qx);
		r2 = this.field.subtract(r1, px);
		r2 = this.field.multiply(s, r2);
		r2 = this.field.add(py, r2);
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
		Polynomial<T> a = r.getEmbedding(this.a);
		Polynomial<T> b = r.getEmbedding(this.b);
		list.add(r.add(r.multiply(3, x, x), r.multiply(a, z, z)));
		list.add(r.multiply(-2, y, z));
		list.add(r.add(r.multiply(2, a, x, z), r.multiply(3, b, z, z), r.multiply(-1, y, y)));
		return list;
	}

	@Override
	public ProjectivePoint<T> getRandomElement() {
		T x;
		T y2;
		do {
			x = this.field.getRandomElement();
			y2 = this.univariateRing.evaluate(this.rhs, Collections.singletonList(x));
		} while (!this.field.hasSqrt(y2));
		T y = this.field.sqrt(y2).keySet().iterator().next();
		return new ProjectivePoint<T>(field, x, y, field.one());

	}

	private ProjectivePoint<T> getRandomElement(int degree) {
		if (!field.isFinite()) {
			throw new RuntimeException("This is wrong!");
		}
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
			T y2 = this.univariateRing.evaluate(this.rhs, Collections.singletonList(x));
			if (!field.hasSqrt(y2)) {
				continue;
			}
			y = field.sqrt(y2).keySet().iterator().next();
			if (y.equals(field.power(y, ch.pow(degree)))) {
				return new ProjectivePoint<T>(field, x, y, field.one());
			}
		}
	}

	@Override
	public boolean isFinite() {
		return this.field.isFinite();
	}

	private void countPointsUpwards(BigInteger basecount, BigInteger q) {
		BigInteger a = BigInteger.TWO;
		BigInteger b = q.add(BigInteger.ONE).subtract(basecount);
		BigInteger base = b;
		int n = 0;
		BigInteger num = field.getNumberOfElements();
		while (!num.equals(BigInteger.ONE)) {
			n++;
			num = num.divide(q);
		}
		for (int i = 1; i < n; i++) {
			BigInteger next = base.multiply(b).subtract(q.multiply(a));
			a = b;
			b = next;
		}
		this.superSingular = true;
		this.numberOfPoints = field.getNumberOfElements().add(BigInteger.ONE).subtract(b);
	}

	private boolean attemptSupersingularCount() {
		this.superSingular = false;
		this.superSingularDeterminied = true;
		T j = this.jInvariant();
		T jp = field.power(j, field.characteristic());
		T jpp = field.power(jp, field.characteristic());
		if (!j.equals(jpp)) {
			return false;
		}
		T ap = field.power(a, field.characteristic());
		T bp = field.power(b, field.characteristic());
		BigInteger one = BigInteger.ONE;
		BigInteger pa1 = field.characteristic().add(one);
		BigInteger pm1 = field.characteristic().subtract(one);
		if (a.equals(ap) && b.equals(bp)) {
			// Curve defined over F_p
			for (int i = 0; i < Math.max(1, 20 - field.characteristic().bitLength() / 2); i++) {
				ProjectivePoint<T> point = this.getRandomElement(1);
				if (!multiply(pa1, point).equals(neutral())) {
					return false;
				}
			}
			countPointsUpwards(pa1, field.characteristic());
			return true;
		}

		T app = field.power(ap, field.characteristic());
		T bpp = field.power(bp, field.characteristic());
		if (!a.equals(app) || !b.equals(bpp)) {
			EllipticCurve<T> newCurve = EllipticCurve.fromJInvariant(field, j);
			this.superSingular = newCurve.attemptSupersingularCount();
			return this.superSingular;
		}
		// Curve defined over F_p^2
		int addOneCount = 0;
		int minusOneCount = 0;
		for (int i = 0; i < Math.max(1, 20 - field.characteristic().bitLength()); i++) {
			ProjectivePoint<T> point = this.getRandomElement(2);

			ProjectivePoint<T> multiplied = multiply(pm1, point);
			if (multiplied.equals(neutral())) {
				minusOneCount++;
				continue;
			}
			multiplied = add(multiplied, add(point, point));
			if (multiply(pa1, point).equals(neutral())) {
				addOneCount++;
				continue;
			}
			return false;
		}
		if (addOneCount > minusOneCount) {
			countPointsUpwards(pa1.multiply(pa1), field.characteristic().pow(2));
		} else {
			countPointsUpwards(pm1.multiply(pm1), field.characteristic().pow(2));
		}
		return true;
	}

	@Override
	public BigInteger getNumberOfElements() {
		if (!this.field.isFinite())
			return BigInteger.valueOf(-1);
		if (this.numberOfPoints.compareTo(BigInteger.ZERO) > 0) {
			return this.numberOfPoints;
		}
		if (field.characteristic().compareTo(BigInteger.valueOf(3)) > 0 && this.attemptSupersingularCount()) {
			return this.getNumberOfElements();
		}
		BigInteger zero = BigInteger.ZERO;
		BigInteger one = BigInteger.ONE;
		BigInteger two = BigInteger.TWO;
		List<BigInteger> primes = new ArrayList<BigInteger>();
		List<BigInteger> t = new ArrayList<BigInteger>();
		BigInteger i = two;
		BigInteger product = one;
		BigInteger q = this.field.getNumberOfElements();
		BigInteger twoSqrtQ = q.sqrt().add(one).shiftLeft(1);
		while (product.compareTo(twoSqrtQ.shiftLeft(1)) < 0) {
			if (i.mod(this.field.characteristic()).equals(zero)) {
				i = i.nextProbablePrime();
				continue;
			}
			if (i.subtract(one).mod(this.field.characteristic()).equals(zero)) {
				i = i.nextProbablePrime();
				continue;
			}
			if (i.add(one).mod(this.field.characteristic()).equals(zero)) {
				i = i.nextProbablePrime();
				continue;
			}
			primes.add(i);
			product = product.multiply(i);
			i = i.nextProbablePrime();
		}
		PolynomialRing<T> r = this.affineRing;
		Polynomial<T> x = r.getVar(1);
		Polynomial<T> y = r.getVar(2);
		PolynomialRing<T> univarring = this.univariateRing;

		for (BigInteger l : primes) {
			if (l.equals(two)) {
				// Computing whether curve has two torsion points.
				if (field.hasRoots(univarring.add(univarring.getVarPower(1, 3),
						univarring.multiply(this.a, univarring.getVar(1)), univarring.getEmbedding(b)))) {
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
			CoordinateRing<T> cr = new CoordinateRing<T>(r, ideal);
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
		CoordinateRingElement<T> xr = cr.multiply(s, s);
		xr = cr.subtract(xr, cr.add(x1, x2));
		CoordinateRingElement<T> yr = cr.add(y1, cr.multiply(s, cr.subtract(xr, x1)));
		List<CoordinateRingElement<T>> result = new ArrayList<>();
		result.add(xr);
		result.add(cr.negative(yr));
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

	public ProjectiveMorphism<T> translationMorphism(ProjectivePoint<T> point) {
		if (point.equals(pointAtInfinity)) {
			List<Polynomial<T>> values = new ArrayList<>();
			values.add(projectiveRing.getVar(1));
			values.add(projectiveRing.getVar(2));
			values.add(projectiveRing.getVar(3));
			return new ProjectiveMorphism<>(asProjectiveVariety(), asProjectiveVariety(), values);
		}
		PolynomialRing<T> r = projectiveRing;
		Polynomial<T> xp = r.getEmbedding(point.getDehomogenisedCoord(1, 3));
		Polynomial<T> yp = r.getEmbedding(point.getDehomogenisedCoord(2, 3));
		Polynomial<T> xmxp = r.subtract(r.getVar(1), r.multiply(xp, r.getVar(3)));
		Polynomial<T> ymyp = r.subtract(r.getVar(2), r.multiply(yp, r.getVar(3)));
		Polynomial<T> x2 = r.getVarPower(1, 2);
		Polynomial<T> xz = r.multiply(r.getVar(1), r.getVar(3));
		Polynomial<T> yz = r.multiply(r.getVar(2), r.getVar(3));
		Polynomial<T> z2 = r.getVarPower(3, 2);
		Polynomial<T> xpx2 = r.multiply(xp, x2);
		Polynomial<T> xp2axz = r.multiply(r.add(r.power(xp, 2), r.getEmbedding(a)), xz);
		Polynomial<T> yyz = r.multiply(-2, yp, yz);
		Polynomial<T> cz2 = r.multiply(r.add(r.multiply(2, r.getEmbedding(b)), r.multiply(a, xp)), z2);
		Polynomial<T> xR = r.add(xpx2, xp2axz, yyz, cz2);
		Polynomial<T> xResult = r.multiply(xmxp, xR);
		Polynomial<T> yResult = r.negative(r.add(r.multiply(yp, r.power(xmxp, 3)),
				r.multiply(ymyp, r.subtract(xR, r.multiply(xp, r.power(xmxp, 2))))));
		Polynomial<T> zResult = r.power(xmxp, 3);
		List<Polynomial<T>> values = new ArrayList<>();
		values.add(xResult);
		values.add(yResult);
		values.add(zResult);
		return new ProjectiveMorphism<>(asProjectiveVariety(), asProjectiveVariety(), values);
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
			return new ProjectiveMorphism<>(asProjectiveVariety(), asProjectiveVariety(), asPolynomials);
		}
		if (n == 1) {
			return new Isomorphism<>(this, negative ? field.getInteger(-1) : field.getInteger(1)).asMorphism();
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
		return new ProjectiveMorphism<>(asProjectiveVariety(), asProjectiveVariety(), polynomials);
	}

	public List<ProjectivePoint<T>> getTorsionPoints(int l) {
		if (l == 1) {
			return Collections.singletonList(this.pointAtInfinity);
		}
		List<ProjectivePoint<T>> torsionPoints = new ArrayList<>(l * l);
		if (l == 2) {
			torsionPoints.add(this.pointAtInfinity);
			Map<T, Integer> rhsRoots = this.field.roots(rhs);
			for (T root : rhsRoots.keySet()) {
				torsionPoints.add(new ProjectivePoint<T>(this.field, root, this.field.zero(), this.field.one()));
			}
		}
		if (l > 2) {
			Polynomial<T> divPoly = null;
			if (l % 2 == 0) {
				int n = l / 2;
				torsionPoints.addAll(this.getTorsionPoints(n));
				for (int i = n - 2; i <= n + 2; i++) {
					this.getDivisionPolynomial(i);
				}
				PolynomialRing<T> r = this.affineRing;
				Polynomial<T> psiNm2 = this.divisionPolynomials.get(n - 2);
				Polynomial<T> psiNm1 = this.divisionPolynomials.get(n - 1);
				Polynomial<T> psiN1 = this.divisionPolynomials.get(n + 1);
				Polynomial<T> psiN2 = this.divisionPolynomials.get(n + 2);
				Polynomial<T> p = r.subtract(r.multiply(psiN2, r.power(psiNm1, 2)),
						r.multiply(psiNm2, r.power(psiN1, 2)));
				divPoly = this.univariateRing.getEmbedding(p, new int[] { 0 });
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

			@SuppressWarnings("unchecked")
			@Override
			public ProjectivePoint<T> next() {
				ProjectivePoint<T> point = nextPoint;
				if (nextNextPoint != null) {
					nextPoint = nextNextPoint;
					nextNextPoint = null;
					return point;
				}
				T x;
				Map<T, Integer> ys;
				do {
					if (!xIterator.hasNext()) {
						nextPoint = null;
						return point;
					}
					x = xIterator.next();
					ys = field.sqrt(univariateRing.evaluate(rhs, x));
				} while (ys.isEmpty());
				Iterator<T> yIterator = ys.keySet().iterator();
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
		String aString = this.a.toString();
		String bString = this.b.toString();
		if (aString.contains(" ")) {
			aString = "(" + aString + ")";
		}
		if (bString.contains(" ")) {
			bString = "(" + bString + ")";
		}
		return "Y^2 = X^3 + " + aString + "X + " + bString;
	}

	@Override
	public List<RationalFunction<T>> getRiemannRochSpace(Divisor<T> div) {
		if (div.getDegree() < 0)
			return Collections.emptyList();
		List<RationalFunction<T>> functions = new ArrayList<>();
		FunctionField<T> ff = getFunctionField();
		List<ProjectivePoint<T>> zeroes = div.getPoles();
		List<ProjectivePoint<T>> poles = div.getZeroes();
		if (poles.size() == 0) {
			return Collections.singletonList(ff.one());}
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
			f = ff.multiply(f, this.getRationalFunction( zero, firstpole, pole));
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
			functions.add(ff.multiply(f, this.getRationalFunction( zero, firstpole, pole)));
		}
		return functions;
	}

	@Override
	public boolean isPrincipal(Divisor<T> div) {
		if (div.getDegree() != 0)
			return false;
		return this.getRiemannRochSpace(div).size() == 1;
	}
}
