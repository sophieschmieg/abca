package fields;

import java.util.List;

import varieties.ProjectivePoint;
import varieties.curves.SmoothCurve;

public class RationalFunction<T extends Element> implements Element {
	private Field<T> field;
	private Polynomial<T> numerator;
	private Polynomial<T> denominator;
	private PolynomialRing<T> ring;
	private SmoothCurve<T> curve;
	
	public static<T extends Element> RationalFunction<T> getEmbedding(T t, SmoothCurve<T> curve, PolynomialRing<T> ring) {
		return new RationalFunction<T>(curve.getField(), ring.getEmbedding(t), ring.getEmbedding(curve.getField().one()), curve, ring);
	}

	public RationalFunction(Field<T> field,
			Polynomial<T> numerator, Polynomial<T> denominator, SmoothCurve<T> curve) {
		this(field, numerator, denominator, curve, new PolynomialRing<T>(field, curve.getEmbeddingDimension() + 1, Polynomial.LEX));
	}
	public RationalFunction(Field<T> field,
			Polynomial<T> numerator, Polynomial<T> denominator, SmoothCurve<T> curve,
			PolynomialRing<T> ring) {
		this.field = field;
		this.ring = ring;
		this.numerator = numerator;
		this.denominator = denominator;
		this.curve = curve;	
		if (this.numerator.getNumVars() != this.denominator.getNumVars())
			throw new ArithmeticException("wrong number of variables");
		if (this.curve.getEmbeddingDimension() + 1 != this.denominator.getNumVars())
			throw new ArithmeticException("wrong number of variables");
		if (this.denominator.equals(this.ring.zero()))
			throw new ArithmeticException("Divison by zero!");
		else if (this.numerator.equals(this.ring.zero()))
			this.denominator = this.ring.one();
		else if (this.numerator.getDegree() != this.denominator.getDegree())
			throw new ArithmeticException("wrong degrees");
		else if (this.denominator.equals(this.numerator)) {
			this.numerator = this.ring.one();
			this.denominator = this.ring.one();
		}
	}
	public ProjectivePoint<T> evaluate(ProjectivePoint<T> p) {
		if (!this.curve.hasRationalPoint(p))
			throw new ArithmeticException("Point not on curve");
		return this.evaluate(p, this.numerator, this.denominator);
	}
	@SuppressWarnings("unchecked")
	private ProjectivePoint<T> evaluate(ProjectivePoint<T> p, Polynomial<T> numerator, Polynomial<T> denominator) {
		T evalnum = numerator.evaluate(p.getCoords());
		T evaldenom = denominator.evaluate(p.getCoords());
		if (!this.field.zero().equals(evalnum) || 
				!this.field.zero().equals(evaldenom))
			return new ProjectivePoint<T>(field, evalnum, evaldenom);
		List<Polynomial<T>> cotangent = this.curve.getCotangentSpace(p);
		Polynomial<T> nextnum = this.ring.zero();
		Polynomial<T> nextdenom = this.ring.zero();
		for (int i = 1; i <= this.curve.getEmbeddingDimension() + 1; i++) {
			nextnum = this.ring.add(nextnum, this.ring.multiply(cotangent.get(i - 1), numerator.derivative(i)));
			nextdenom = this.ring.add(nextdenom, this.ring.multiply(cotangent.get(i - 1), denominator.derivative(i)));
		}
		return this.evaluate(p, nextnum, nextdenom);
	}
	public Polynomial<T> getNumerator() {
		return numerator;
	}
	public Polynomial<T> getDenominator() {
		return denominator;
	}
	@Override
	public boolean equals(Object O) {
		if (!(O instanceof RationalFunction))
			return false;
		@SuppressWarnings("unchecked")
		RationalFunction<T> f = (RationalFunction<T>)O;
		return this.ring.multiply(this.numerator, f.denominator).equals(this.ring.multiply(f.numerator, this.denominator));
	}
	@Override
	public String toString() {
		return "(" + this.numerator.toString() + ")/(" + this.denominator.toString() + ")";
	}
	public int getOrder(ProjectivePoint<T> p) {
		return this.getOrder(p, this.numerator) - this.getOrder(p, this.denominator);
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
		RationalFunction<T> f = (RationalFunction<T>)o;
		return this.ring.multiply(this.numerator, f.denominator).compareTo(this.ring.multiply(f.numerator, this.denominator));
	}
}
