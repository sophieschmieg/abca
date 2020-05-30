package varieties.curves;

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;

import fields.PrimeField;
import fields.PrimeField.PrimeFieldElement;
import fields.CoordinateRing;
import fields.CoordinateRing.CoordinateRingElement;
import fields.Element;
import fields.Field;
import fields.FunctionField;
import fields.Group;
import fields.Polynomial;
import fields.PolynomialRing;
import fields.RationalFunction;
import util.MiscAlgorithms;
import varieties.ProjectivePoint;
import varieties.ProjectiveSpace;
import varieties.ProjectiveVariety;
import varieties.curves.DivisorGroup.Divisor;

public class EllipticCurve<T extends Element> extends ProjectiveVariety<T>
		implements SmoothCurve<T>, Group<ProjectivePoint<T>> {
	private Field<T> field;
	private T a;
	private T b;
	private ProjectivePoint<T> pointAtInfinity;
	private PolynomialRing<T> projectiveRing;
	private PolynomialRing<T> affineRing;
	private PolynomialRing<T> univariateRing;
	private CoordinateRing<T> coordinateRing;
	private Polynomial<T> rhs;
	@SuppressWarnings("rawtypes")
	private static Map<Field, ProjectiveSpace> spaces = new HashMap<Field, ProjectiveSpace>();
	private Map<Integer, Polynomial<T>> divisionPolynomials;
	private Polynomial<T> definingPolynomial;
	private BigInteger numberOfPoints = BigInteger.valueOf(-1);

	@SuppressWarnings("unchecked")
	private static <T extends Element> ProjectiveSpace<T> getSpace(Field<T> field) {
		if (!EllipticCurve.spaces.containsKey(field))
			EllipticCurve.spaces.put(field, new ProjectiveSpace<T>(field, 2));
		return EllipticCurve.spaces.get(field);
	}

	private static <T extends Element> PolynomialRing<T>.Ideal getIdeal(Field<T> field, T a, T b) {
		PolynomialRing<T> ring = EllipticCurve.getSpace(field).getRing();
		Polynomial<T> f = ring.getEmbedding(field.one(), new int[] { 3, 0, 0 });
		f = ring.add(f, ring.getEmbedding(field.negative(field.one()), new int[] { 0, 2, 1 }));
		f = ring.add(f, ring.getEmbedding(a, new int[] { 1, 0, 2 }));
		f = ring.add(f, ring.getEmbedding(b, new int[] { 0, 0, 3 }));
		return ring.getIdeal(Collections.singletonList(f));
	}

	public EllipticCurve(Field<T> field, T a, T b) {
		super(EllipticCurve.getSpace(field), EllipticCurve.getIdeal(field, a, b));
		this.field = field;
		this.a = a;
		this.b = b;
		this.pointAtInfinity = new ProjectivePoint<T>(this.field, this.field.zero(), this.field.one(),
				this.field.zero());
		if (this.field
				.add(this.field.multiply(4, this.field.power(a, 3)), this.field.multiply(27, this.field.power(b, 2)))
				.equals(this.field.zero()))
			throw new ArithmeticException("Singular curve");
		this.projectiveRing = this.getSpace().getRing();
		this.affineRing = new PolynomialRing<T>(this.field, 2, Polynomial.REVLEX);
		PolynomialRing<T> r = this.affineRing;
		Polynomial<T> x = r.getVar(1);
		Polynomial<T> x3 = r.getVar(1, 3);
		Polynomial<T> y2 = r.getVar(2, 2);
		Polynomial<T> rb = r.getEmbedding(this.b);
		Polynomial<T> rhs = r.add(x3, r.multiply(this.a, x), rb);
		this.definingPolynomial = r.subtract(rhs, y2);
		this.univariateRing = new PolynomialRing<T>(this.field, 1, Polynomial.LEX);
		r = this.univariateRing;
		x = r.getVar(1);
		x3 = r.getVar(1, 3);
		rb = r.getEmbedding(this.b);
		this.rhs = r.add(x3, r.multiply(this.a, x), rb);
	}

	public T jInvariant() {
		T fourACube = this.field.multiply(4, this.field.power(a, 3));
		T twentySevenBSquare = this.field.multiply(27, this.field.power(b, 2));
		return this.field.multiply(1728, this.field.divide(fourACube, this.field.add(fourACube, twentySevenBSquare)));
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
		return this.definingPolynomial.evaluate(p.getDehomogenous(3).getCoords()).equals(this.field.zero());
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
		Polynomial<T> x2 = r.getVar(1, 2);
		Polynomial<T> x3 = r.getVar(1, 3);
		Polynomial<T> x4 = r.getVar(1, 4);
		Polynomial<T> x6 = r.getVar(1, 6);
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
	public boolean isProjective() {
		return true;
	}

	public RationalFunction<T> getRationalFunction(ProjectivePoint<T> zero, ProjectivePoint<T> pole1,
			ProjectivePoint<T> pole2) {
		ProjectivePoint<T> common = this.getThirdIntersection(pole1, pole2);
		Polynomial<T> numerator;
		Polynomial<T> denominator;
		if (common.equals(zero))
			numerator = this.getTangentSpace(zero).get(0);
		else
			numerator = this.getSpace().asHyperplaneIdeal(zero, common).getBasis().first();
		if (pole1.equals(pole2))
			denominator = this.getTangentSpace(pole1).get(0);
		else
			denominator = this.getSpace().asHyperplaneIdeal(pole1, pole2).getBasis().first();
		return new RationalFunction<T>(this.field, numerator, denominator, this, this.projectiveRing);
	}

	public FunctionField<T> getFunctionField() {
		return new FunctionField<T>(field, this, this.projectiveRing);
	}

	private List<Polynomial<T>> splitCoordinateRingElement(CoordinateRingElement<T> t) {
		List<T> values = new ArrayList<T>();
		values.add(null);
		values.add(field.zero());
		Polynomial<T> xPart = t.getElement().evaluatePartially(values);
		values.clear();
		values.add(null);
		values.add(field.one());
		Polynomial<T> yPart = affineRing.subtract(t.getElement().evaluatePartially(values), xPart);
		List<Polynomial<T>> result = new ArrayList<>();
		result.add(univariateRing.getEmbedding(xPart, new int[] { 0 }));
		result.add(univariateRing.getEmbedding(yPart, new int[] { 0 }));
		return result;
	}

	private CoordinateRingElement<T> joinCoordinateRingElement(Polynomial<T> xPart, Polynomial<T> yPart) {
		CoordinateRing<T> cr = getCoordinateRing();
		PolynomialRing<T> r = affineRing;
		CoordinateRingElement<T> y = cr.getEmbedding(affineRing.getVar(2));
		CoordinateRingElement<T> affineX = cr.getEmbedding(r.getEmbedding(xPart, new int[] { 0 }));
		CoordinateRingElement<T> affineY = cr.getEmbedding(r.getEmbedding(yPart, new int[] { 0 }));
		return cr.add(affineX, cr.multiply(y, affineY));
	}

	@Override
	public boolean hasSimplify() {
		return true;
	}

	@Override
	public List<CoordinateRingElement<T>> simplify(RationalFunction<T> t) {
		CoordinateRing<T> cr = getCoordinateRing();
		PolynomialRing<T> r = univariateRing;
		CoordinateRingElement<T> num = t.getNumerator();
		CoordinateRingElement<T> denom = t.getDenominator();
		List<Polynomial<T>> splitDenom = splitCoordinateRingElement(t.getDenominator());
		CoordinateRingElement<T> normalizer = joinCoordinateRingElement(splitDenom.get(0),
				r.negative(splitDenom.get(1)));
		num = cr.multiply(normalizer, num);
		denom = cr.multiply(normalizer, denom);
		List<Polynomial<T>> splitNum = splitCoordinateRingElement(num);
		Polynomial<T> numX = splitNum.get(0);
		Polynomial<T> numY = splitNum.get(1);
		splitDenom = splitCoordinateRingElement(denom);
		Polynomial<T> denomX = splitDenom.get(0);
		Polynomial<T> denomY = splitDenom.get(1);
		if (!denomY.equals(r.zero())) {
			throw new ArithmeticException("Could not simplify, normalizer failed!");
		}
		Polynomial<T> gcd = r.gcd(numX, numY);
		gcd = r.gcd(gcd, denomX);
		numX = r.quotientAndRemainder(numX, gcd).get(0);
		numY = r.quotientAndRemainder(numY, gcd).get(0);
		denomX = r.quotientAndRemainder(denomX, gcd).get(0);
		List<CoordinateRingElement<T>> result = new ArrayList<>();
		result.add(joinCoordinateRingElement(numX, numY));
		result.add(joinCoordinateRingElement(denomX, r.zero()));
		return result;
	}

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
			line.add(poly.evaluate(p.getCoords()));
		return Collections.singletonList(this.projectiveRing.getLinear(line));
	}

	@Override
	public List<Polynomial<T>> getCotangentSpace(ProjectivePoint<T> p) {
		List<Polynomial<T>> list = this.getDifferentials();
		List<Polynomial<T>> reslist = new ArrayList<Polynomial<T>>();
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
		List<Polynomial<T>> list = new ArrayList<Polynomial<T>>();
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
			y2 = this.rhs.evaluate(x);
		} while (!this.field.hasSqrt(y2));
		T y = this.field.sqrt(y2).iterator().next();
		return new ProjectivePoint<T>(field, x, y, field.one());

	}

	private ProjectivePoint<T> getRandomElement(int degree) {
		if (!field.isFinite()) {
			throw new RuntimeException("This is wrong!");
		}
		T x;
		T y2;
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
		do {
			T xq = this.field.getRandomElement();
			x = field.zero();
			for (int i = 0; i < n; i += degree) {
				x = field.add(x, field.power(xq, ch.pow(i)));
			}
			y2 = this.rhs.evaluate(x);
		} while (!this.field.hasSqrt(y2));
		T y = this.field.sqrt(y2).iterator().next();
		return new ProjectivePoint<T>(field, x, y, field.one());

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
		this.numberOfPoints = field.getNumberOfElements().add(BigInteger.ONE).subtract(b);
	}

	private boolean attemptSupersingularCount() {
		T j = this.jInvariant();
		T jp = field.power(j, field.characteristic());
		T jpp = field.power(jp, field.characteristic());
		if (!j.equals(jpp)) {
			return false;
		}
		BigInteger one = BigInteger.ONE;
		BigInteger pa1 = field.characteristic().add(one);
		if (j.equals(jp)) {
			ProjectivePoint<T> p = this.getRandomElement(1);
			if (multiply(pa1, p).equals(neutral())) {
				countPointsUpwards(pa1, field.characteristic());
				return true;
			} else {
				return false;
			}
		} else {
			BigInteger pm1 = field.characteristic().subtract(one);
			ProjectivePoint<T> p = this.getRandomElement(2);
			if (multiply(pa1, p).equals(neutral())) {
				countPointsUpwards(pa1.pow(2), field.characteristic().pow(2));
				return true;
			} else if (multiply(pm1, p).equals(neutral())) {
				countPointsUpwards(pm1.pow(2), field.characteristic().pow(2));
				return true;
			} else {
				return false;
			}
		}
	}

	@Override
	public BigInteger getNumberOfElements() {
		if (!this.field.isFinite())
			return BigInteger.valueOf(-1);
		if (this.numberOfPoints.compareTo(BigInteger.ZERO) > 0) {
			return this.numberOfPoints;
		}
		if (field.characteristic().compareTo(BigInteger.valueOf(100)) > 0 && this.attemptSupersingularCount()) {
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
				if (univarring.hasRoots(univarring.add(univarring.getVar(1, 3),
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
			idealGen.add(psiL);
			idealGen.add(this.definingPolynomial);
			PolynomialRing<T>.Ideal ideal = r.getIdeal(idealGen);
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
				PrimeField lField = new PrimeField(l);
				PrimeFieldElement qModL = lField.getInteger(q);
				if (lField.hasSqrt(qModL)) {
					PrimeFieldElement w = lField.sqrt(qModL).iterator().next();
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
		List<Polynomial<T>> egcd = univarring.extendedEuclidean(univarA, psiL);
		Polynomial<T> univarAInv = univarring.multiply(egcd.get(1), univarring.inverse(egcd.get(0)));
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
		xdiffPoly = r.quotientAndRemainder(xdiffPoly, gcdPoly).get(0);
		List<Polynomial<T>> ydiffqr = this.affineRing.quotientAndRemainder(ydiff.getElement(), gcdMulti);
		if (!ydiffqr.get(1).equals(this.affineRing.zero())) {
			throw new RuntimeException("Cannot add generic!");
		}
		CoordinateRingElement<T> ydiffReduced = cr.getEmbedding(ydiffqr.get(0));
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
		psiNsq = r.reduce(psiNsq, Collections.singletonList(this.definingPolynomial));
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

	public List<ProjectivePoint<T>> getTorsionPoints(int l) {
		if (l == 1) {
			return Collections.singletonList(this.pointAtInfinity);
		}
		List<ProjectivePoint<T>> torsionPoints = new ArrayList<>(l * l);
		if (l == 2) {
			torsionPoints.add(this.pointAtInfinity);
			List<T> rhsRoots = this.univariateRing.roots(rhs);
			for (T root : rhsRoots) {
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
			List<T> xs = this.univariateRing.roots(divPoly);
			for (T x : xs) {
				for (T y : this.field.sqrt(this.rhs.evaluate(Collections.singletonList(x)))) {
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
				Set<T> ys;
				do {
					if (!xIterator.hasNext()) {
						nextPoint = null;
						return point;
					}
					x = xIterator.next();
					ys = field.sqrt(rhs.evaluate(Collections.singletonList(x)));
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

	@Override
	public String toString() {
		return "Y^2Z = X^3 + " + this.a.toString() + "X + " + this.b.toString();
	}

	@Override
	public List<RationalFunction<T>> getRiemannRochSpace(Divisor<T> div) {
		if (div.getDegree() < 0)
			return Collections.emptyList();
		List<RationalFunction<T>> functions = new ArrayList<RationalFunction<T>>();
		FunctionField<T> ff = new FunctionField<T>(this.field, this, this.projectiveRing);
		List<ProjectivePoint<T>> zeroes = div.getPoles();
		List<ProjectivePoint<T>> poles = div.getZeroes();
		if (poles.size() == 0)
			return Collections.singletonList(ff.one());
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
	public boolean isPrincipal(Divisor<T> div) {
		if (div.getDegree() != 0)
			return false;
		return this.getRiemannRochSpace(div).size() == 1;
	}
}
