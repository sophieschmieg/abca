package varieties;

import java.util.ArrayList;
import java.util.List;

import fields.helper.CoordinateRing;
import fields.helper.CoordinateRing.CoordinateRingElement;
import fields.interfaces.Element;
import fields.interfaces.Field;
import fields.interfaces.Ideal;
import fields.interfaces.Polynomial;
import fields.interfaces.PolynomialRing;
import fields.polynomials.AbstractPolynomialRing;
import fields.polynomials.Monomial;
import varieties.curves.SmoothCurve;

public class RationalFunction<T extends Element<T>> implements Element<RationalFunction<T>> {
	private Field<T> field;
	private int numVars;
	private CoordinateRingElement<T> numerator;
	private CoordinateRingElement<T> denominator;
	private PolynomialRing<T> ring;
	private SmoothCurve<T> curve;
	private CoordinateRing<T> coordRing;

	public static <T extends Element<T>> RationalFunction<T> getEmbedding(T t, SmoothCurve<T> curve,
			PolynomialRing<T> ring) {
		return new RationalFunction<T>(curve.getField(), ring.getEmbedding(t),
				ring.getEmbedding(curve.getField().one()), curve, ring);
	}

	public RationalFunction(Field<T> field, Polynomial<T> numerator, Polynomial<T> denominator, SmoothCurve<T> curve) {
		this(field, numerator, denominator, curve,
				AbstractPolynomialRing.getPolynomialRing(field, curve.getEmbeddingDimension() + 1, Monomial.GREVLEX));
	}

	public RationalFunction(Field<T> field, Polynomial<T> numerator, Polynomial<T> denominator, SmoothCurve<T> curve,
			PolynomialRing<T> ring) {
		this.field = field;
		this.ring = ring;
		this.curve = curve;
		this.coordRing = curve.getCoordinateRing();
		PolynomialRing<T> affineRing = coordRing.getPolynomialRing();
		this.numVars = affineRing.numberOfVariables();
		if (numerator.numberOfVariables() == this.curve.getEmbeddingDimension() + 1) {
			if (denominator.numberOfVariables() == this.curve.getEmbeddingDimension() + 1) {
				if (numerator.degree() != denominator.degree() || !ring.isHomogeneous(numerator)
						|| !ring.isHomogeneous(denominator)) {
					throw new ArithmeticException("wrong degrees: " + numerator + "/" + denominator);
				}
				int[] map = new int[numVars];
				for (int i = 0; i < map.length; i++) {
					map[i] = i;
				}
				this.numerator = coordRing
						.getEmbedding(affineRing.getEmbedding(ring.dehomogenize(numerator, numVars + 1), map));
				this.denominator = coordRing
						.getEmbedding(affineRing.getEmbedding(ring.dehomogenize(denominator, numVars + 1), map));
			} else {
				throw new ArithmeticException("wrong number of variables");
			}
		} else if (numerator.numberOfVariables() == this.curve.getEmbeddingDimension()) {
			if (numerator.numberOfVariables() == this.curve.getEmbeddingDimension()) {
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
		this.numVars = this.coordRing.getPolynomialRing().numberOfVariables();
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
		Ideal<Polynomial<T>> gcdIdeal = r.getIdeal(gcdGen);
		if (gcdIdeal.generators().size() == 1) {
			Polynomial<T> gcd = gcdIdeal.generators().get(0);
			numerator = cr.getEmbedding(r.quotientAndRemainder(numerator.getElement(), gcd).get(0));
			denominator = cr.getEmbedding(r.quotientAndRemainder(denominator.getElement(), gcd).get(0));
		}
		List<Polynomial<T>> basis = cr.getIdeal().generators();
		List<Polynomial<T>> gcdIdealGen = new ArrayList<Polynomial<T>>();
		gcdIdealGen.addAll(basis);
		gcdIdealGen.add(numerator.getElement());
		gcdIdealGen.add(denominator.getElement());
		Ideal<Polynomial<T>> gcdBasisIdeal = r.getIdeal(gcdIdealGen);
		List<Polynomial<T>> gcdBasis = new ArrayList<>(gcdBasisIdeal.generators());
		CoordinateRingElement<T> reducedNum = numerator;
		CoordinateRingElement<T> reducedDenom = denominator;
		if (gcdBasis.size() == 2) {
			basisLoop: for (Polynomial<T> b : cr.getIdeal().generators()) {
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
		T norm = reducedDenom.getElement().leadingCoefficient();
		this.numerator = cr.multiply(field.inverse(norm), reducedNum);
		this.denominator = cr.multiply(field.inverse(norm), reducedDenom);
	}

	private List<CoordinateRingElement<T>> simplify(CoordinateRing<T> cr) {
		PolynomialRing<T> r = cr.getPolynomialRing();
		PolynomialRing<T> elimR = r.addVariableWithElimination(1);
		List<Polynomial<T>> gen = new ArrayList<>();
		Polynomial<T> t = elimR.getVar(1);
		Polynomial<T> num = elimR.getEmbeddingWithElimination(numerator.getElement(), 1);
		Polynomial<T> denom = elimR.getEmbeddingWithElimination(denominator.getElement(), 1);
		Polynomial<T> fun = elimR.subtract(num, elimR.multiply(t, denom));
		gen.add(fun);
		for (Polynomial<T> b : cr.getIdeal().generators()) {
			gen.add(elimR.getEmbeddingWithElimination(b, 1));
		}
		Ideal<Polynomial<T>> elimIdeal = elimR.getIdeal(gen);
		Polynomial<T> eliminated = elimIdeal.generators().get(0);
		List<T> eval = new ArrayList<T>();
		eval.add(field.zero());
		for (int i = 0; i < r.numberOfVariables(); i++) {
			eval.add(null);
		}
		Polynomial<T> newNum = r.negative(r.getEmbeddingWithElimination(elimR.partiallyEvaluate(eliminated, eval), -1));
		eval.set(0, field.one());
		Polynomial<T> newDenom = r.add(r.getEmbeddingWithElimination(elimR.partiallyEvaluate(eliminated, eval), -1),
				newNum);
		List<CoordinateRingElement<T>> result = new ArrayList<>();
		result.add(cr.getEmbedding(newNum));
		result.add(cr.getEmbedding(newDenom));
		return result;
	}

	private List<Polynomial<T>> homogenize() {
		Polynomial<T> numerator = this.numerator.getElement();
		Polynomial<T> denominator = this.denominator.getElement();
		int degree = Math.max(numerator.degree(), denominator.degree());
		numerator = ring.multiply(ring.homogenize(numerator),
				ring.getVarPower(numVars + 1, degree - numerator.degree()));
		denominator = ring.multiply(ring.homogenize(denominator),
				ring.getVarPower(numVars + 1, degree - denominator.degree()));
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
			num = r.getEmbedding(ring.dehomogenize(r.homogenize(num), dehom), map);
			denom = r.getEmbedding(ring.dehomogenize(r.homogenize(denom), dehom), map);
			List<Polynomial<T>> gen = new ArrayList<>();
			for (Polynomial<T> b : coordRing.getIdeal().generators()) {
				gen.add(r.getEmbedding(ring.dehomogenize(r.homogenize(b), dehom), map));
			}
			coordR = new CoordinateRing<>(r, r.getIdeal(gen));
		}

		Ideal<Polynomial<T>> local = p.getDehomogenous(dehom).asIdeal(r);
		Ideal<Polynomial<T>> pointIdeal = r.add(coordR.getIdeal(), local);
		CoordinateRing<T> pointRing = new CoordinateRing<>(coordR.getPolynomialRing(), pointIdeal);
		List<CoordinateRingElement<T>> simplified = simplify(pointRing);
		if (simplified.get(0).getElement().degree() > 0 || simplified.get(1).getElement().degree() > 0) {
			throw new ArithmeticException("Something went wrong!");
		}
		// System.err.println("simp " + simplified);
		T first = simplified.get(0).getElement().leadingCoefficient();
		T second = simplified.get(1).getElement().leadingCoefficient();
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
		T eval = this.ring.evaluate(poly, p.getCoords());
		if (!this.field.zero().equals(eval))
			return 0;
		List<Polynomial<T>> cotangent = this.curve.getCotangentSpace(p);
		Polynomial<T> next = this.ring.zero();
		for (int i = 1; i <= this.curve.getEmbeddingDimension() + 1; i++) {
			next = this.ring.add(next, this.ring.multiply(cotangent.get(i - 1), this.ring.derivative(poly, i)));
		}
		return this.getOrder(p, next) + 1;
	}

	@Override
	public int compareTo(RationalFunction<T> f) {
		return this.coordRing.multiply(this.numerator, f.denominator)
				.compareTo(this.coordRing.multiply(f.numerator, this.denominator));
	}
}
