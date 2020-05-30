package fields;

import java.util.ArrayList;
import java.util.List;
import java.util.SortedSet;

import fields.CoordinateRing.CoordinateRingElement;
import varieties.ProjectivePoint;
import varieties.curves.SmoothCurve;

public class RationalFunction<T extends Element> implements Element {
	private Field<T> field;
	private int numVars;
	private CoordinateRingElement<T> numerator;
	private CoordinateRingElement<T> denominator;
	private PolynomialRing<T> ring;
	private SmoothCurve<T> curve;
	private CoordinateRing<T> coordRing;

	public static <T extends Element> RationalFunction<T> getEmbedding(T t, SmoothCurve<T> curve,
			PolynomialRing<T> ring) {
		return new RationalFunction<T>(curve.getField(), ring.getEmbedding(t),
				ring.getEmbedding(curve.getField().one()), curve, ring);
	}

	public RationalFunction(Field<T> field, Polynomial<T> numerator, Polynomial<T> denominator, SmoothCurve<T> curve) {
		this(field, numerator, denominator, curve,
				new PolynomialRing<T>(field, curve.getEmbeddingDimension() + 1, Polynomial.GREVLEX));
	}

	public RationalFunction(Field<T> field, Polynomial<T> numerator, Polynomial<T> denominator, SmoothCurve<T> curve,
			PolynomialRing<T> ring) {
		this.field = field;
		this.ring = ring;
		this.curve = curve;
		this.coordRing = curve.getCoordinateRing();
		PolynomialRing<T> affineRing = coordRing.getPolynomialRing();
		this.numVars = affineRing.getNumVars();
		if (numerator.getNumVars() == this.curve.getEmbeddingDimension() + 1) {
			if (denominator.getNumVars() == this.curve.getEmbeddingDimension() + 1) {
				if (numerator.getDegree() != denominator.getDegree() || !numerator.isHomogeneous()
						|| !denominator.isHomogeneous()) {
					throw new ArithmeticException("wrong degrees: " + numerator + "/" + denominator);
				}
				int[] map = new int[numVars];
				for (int i = 0; i < map.length; i++) {
					map[i] = i;
				}
				this.numerator = coordRing.getEmbedding(affineRing.getEmbedding(numerator.dehom(numVars + 1), map));
				this.denominator = coordRing.getEmbedding(affineRing.getEmbedding(denominator.dehom(numVars + 1), map));
			} else {
				throw new ArithmeticException("wrong number of variables");
			}
		} else if (numerator.getNumVars() == this.curve.getEmbeddingDimension()) {
			if (numerator.getNumVars() == this.curve.getEmbeddingDimension()) {
				this.numerator = coordRing.getEmbedding(numerator);
				this.denominator = coordRing.getEmbedding(denominator);

			} else {
				throw new ArithmeticException("wrong number of variables");
			}
		} else {
			throw new ArithmeticException("wrong number of variables");
		}
	}

	public RationalFunction(Field<T> field, CoordinateRingElement<T> numerator, CoordinateRingElement<T> denominator,
			SmoothCurve<T> curve, PolynomialRing<T> ring) {
		this.field = field;
		this.curve = curve;
		this.ring = ring;
		this.coordRing = curve.getCoordinateRing();
		this.numVars = this.coordRing.getPolynomialRing().getNumVars();
		this.numerator = numerator;
		this.denominator = denominator;
	}

	public void simplify() {
		if (curve.hasSimplify()) {
			List<CoordinateRingElement<T>> simplified = curve.simplify(this);
			this.numerator = simplified.get(0);
			this.denominator = simplified.get(1);
			return;
		}
		CoordinateRing<T> cr = coordRing;
		List<Polynomial<T>> gcdGen = new ArrayList<>();
		gcdGen.add(numerator.getElement());
		gcdGen.add(denominator.getElement());
		PolynomialRing<T> r = cr.getPolynomialRing();
		PolynomialRing<T>.Ideal gcdIdeal = r.getIdeal(gcdGen);
		if (gcdIdeal.getBasis().size() == 1) {
			Polynomial<T> gcd = gcdIdeal.getBasis().first();
			numerator = cr.getEmbedding(r.quotientAndRemainder(numerator.getElement(), gcd).get(0));
			denominator = cr.getEmbedding(r.quotientAndRemainder(denominator.getElement(), gcd).get(0));
		}
		SortedSet<Polynomial<T>> basis = cr.getIdeal().getBasis();
		List<Polynomial<T>> gcdIdealGen = new ArrayList<Polynomial<T>>();
		gcdIdealGen.addAll(basis);
		gcdIdealGen.add(numerator.getElement());
		gcdIdealGen.add(denominator.getElement());
		PolynomialRing<T>.Ideal gcdBasisIdeal = r.getIdeal(gcdIdealGen);
		List<Polynomial<T>> gcdBasis = new ArrayList<>(gcdBasisIdeal.getBasis());
		CoordinateRingElement<T> reducedNum = numerator;
		CoordinateRingElement<T> reducedDenom = denominator;
		if (gcdBasis.size() == 2) {
			basisLoop: for (Polynomial<T> b : cr.getIdeal().getBasis()) {
				List<Polynomial<T>> reducedBasis = r.generalQuotientAndRemainder(b, gcdBasis).subList(0, 2);
				System.err.println("reduced Basis" + b + " = " + reducedBasis);
				System.err.println("gcd Basis" + gcdBasis);
				for (int i = 0; i < 2; i++) {
					Polynomial<T> coeff = reducedBasis.get(i);
					if (r.isUnit(coeff)) {
						Polynomial<T> other = r.multiply(r.inverse(r.negative(coeff)), reducedBasis.get(1 - i));
						List<Polynomial<T>> numCoeffs = r.generalQuotientAndRemainder(reducedNum.getElement(),
								gcdBasis);
						System.err.println(numCoeffs);
						reducedNum = cr.getEmbedding(numCoeffs.get(1 - i));
						reducedNum = cr.add(reducedNum, cr.getEmbedding(r.multiply(other, numCoeffs.get(i))));
						List<Polynomial<T>> denomCoeffs = r.generalQuotientAndRemainder(reducedDenom.getElement(),
								gcdBasis);
						reducedDenom = cr.getEmbedding(denomCoeffs.get(1 - i));
						reducedDenom = cr.add(reducedDenom, cr.getEmbedding(r.multiply(other, denomCoeffs.get(i))));
						break basisLoop;
					}
				}
			}
		}
		T norm = reducedDenom.getElement().getLeadingCoefficient();
		this.numerator = cr.multiply(field.inverse(norm), reducedNum);
		this.denominator = cr.multiply(field.inverse(norm), reducedDenom);
	}

	private List<CoordinateRingElement<T>> simplify(CoordinateRing<T> cr) {
		PolynomialRing<T> r = cr.getPolynomialRing();
		PolynomialRing<T> elimR = r.addVariableWithElimination(1);
		List<Polynomial<T>> gen = new ArrayList<>();
		Polynomial<T> t = elimR.getVar(1);
		Polynomial<T> num = elimR.getEmbeddingShift(numerator.getElement(), 1);
		Polynomial<T> denom = elimR.getEmbeddingShift(denominator.getElement(), 1);
		Polynomial<T> fun = elimR.subtract(num, elimR.multiply(t, denom));
		gen.add(fun);
		for (Polynomial<T> b : cr.getIdeal().getBasis()) {
			gen.add(elimR.getEmbeddingShift(b, 1));
		}
		PolynomialRing<T>.Ideal elimIdeal = elimR.getIdeal(gen);
		Polynomial<T> eliminated = elimIdeal.getBasis().first();
		List<T> eval = new ArrayList<T>();
		eval.add(field.zero());
		for (int i = 0; i < r.getNumVars(); i++) {
			eval.add(null);
		}
		Polynomial<T> newNum = r.negative(r.getEmbeddingShift(eliminated.evaluatePartially(eval), -1));
		eval.set(0, field.one());
		Polynomial<T> newDenom = r.add(r.getEmbeddingShift(eliminated.evaluatePartially(eval), -1), newNum);
		List<CoordinateRingElement<T>> result = new ArrayList<>();
		result.add(cr.getEmbedding(newNum));
		result.add(cr.getEmbedding(newDenom));
		return result;
	}

	private List<Polynomial<T>> homogenize() {
		Polynomial<T> numerator = this.numerator.getElement();
		Polynomial<T> denominator = this.denominator.getElement();
		int degree = Math.max(numerator.getDegree(), denominator.getDegree());
		numerator = ring.multiply(numerator.homogenize(), ring.getVar(numVars + 1, degree - numerator.getDegree()));
		denominator = ring.multiply(denominator.homogenize(),
				ring.getVar(numVars + 1, degree - denominator.getDegree()));
		List<Polynomial<T>> result = new ArrayList<>();
		result.add(numerator);
		result.add(denominator);
		return result;
	}

	public ProjectivePoint<T> evaluate(ProjectivePoint<T> p) {
		if (!this.curve.hasRationalPoint(p))
			throw new ArithmeticException("Point not on curve");
		int dehom = numVars + 1;
		CoordinateRing<T> coordR = coordRing;
		PolynomialRing<T> r = coordR.getPolynomialRing();
		Polynomial<T> num = numerator.getElement();
		Polynomial<T> denom = denominator.getElement();
		if (p.getCoord(dehom).equals(field.zero())) {
			dehom = p.getNonZero();
			int[] map = new int[numVars + 1];
			for (int i = 0; i < dehom - 1; i++) {
				map[i] = i;
			}
			map[dehom - 1] = -1;
			for (int i = dehom - 1; i < numVars; i++) {
				map[i + 1] = i;
			}
			num = r.getEmbedding(num.homogenize().dehom(dehom), map);
			denom = r.getEmbedding(denom.homogenize().dehom(dehom), map);
			List<Polynomial<T>> gen = new ArrayList<>();
			for (Polynomial<T> b : coordRing.getIdeal().getBasis()) {
				gen.add(r.getEmbedding(b.homogenize().dehom(dehom), map));
			}
			coordR = new CoordinateRing<>(r, r.getIdeal(gen));
		}

		PolynomialRing<T>.Ideal local = p.getDehomogenous(dehom).asIdeal(r);
		PolynomialRing<T>.Ideal pointIdeal = coordR.getIdeal().add(local);
		CoordinateRing<T> pointRing = new CoordinateRing<T>(coordR.getPolynomialRing(), pointIdeal);
		List<CoordinateRingElement<T>> simplified = simplify(pointRing);
		if (simplified.get(0).getElement().getDegree() > 0 || simplified.get(1).getElement().getDegree() > 0) {
			throw new ArithmeticException("Something went wrong!");
		}
		// System.err.println("simp " + simplified);
		T first = simplified.get(0).getElement().getLeadingCoefficient();
		T second = simplified.get(1).getElement().getLeadingCoefficient();
		if (first.equals(field.zero()) && second.equals(field.zero())) {
			RationalFunction<T> f = new RationalFunction<>(field, denominator, numerator, curve, ring);
			f.simplify();
			ProjectivePoint<T> result = f.evaluate(p);
			return new ProjectivePoint<T>(field, result.getCoord(2), result.getCoord(1));
		}
		return new ProjectivePoint<T>(field, first, second);
	}

	public CoordinateRingElement<T> getNumerator() {
		return numerator;
	}

	public CoordinateRingElement<T> getDenominator() {
		return denominator;
	}

	@Override
	public boolean equals(Object O) {
		if (!(O instanceof RationalFunction))
			return false;
		@SuppressWarnings("unchecked")
		RationalFunction<T> f = (RationalFunction<T>) O;
		return this.coordRing.multiply(this.numerator, f.denominator)
				.equals(this.coordRing.multiply(f.numerator, this.denominator));
	}

	@Override
	public String toString() {
		return "(" + this.numerator.toString() + ")/(" + this.denominator.toString() + ")";
	}

	public int getOrder(ProjectivePoint<T> p) {
		List<Polynomial<T>> homogenous = this.homogenize();
		return this.getOrder(p, homogenous.get(0)) - this.getOrder(p, homogenous.get(1));
	}

	private int getOrder(ProjectivePoint<T> p, Polynomial<T> poly) {
		T eval = poly.evaluate(p.getCoords());
		if (!this.field.zero().equals(eval))
			return 0;
		List<Polynomial<T>> cotangent = this.curve.getCotangentSpace(p);
		Polynomial<T> next = this.ring.zero();
		for (int i = 1; i <= this.curve.getEmbeddingDimension() + 1; i++) {
			next = this.ring.add(next, this.ring.multiply(cotangent.get(i - 1), poly.derivative(i)));
		}
		return this.getOrder(p, next) + 1;
	}

	@Override
	public int compareTo(Element o) {
		@SuppressWarnings("unchecked")
		RationalFunction<T> f = (RationalFunction<T>) o;
		return this.coordRing.multiply(this.numerator, f.denominator)
				.compareTo(this.coordRing.multiply(f.numerator, this.denominator));
	}
}
