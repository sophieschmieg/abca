package fields;

import varieties.curves.SmoothCurve;


public class FunctionField<T extends Element> extends AbstractField<RationalFunction<T>> {
	private Field<T> field;
	private PolynomialRing<T> ring;
	private SmoothCurve<T> curve;
	public FunctionField(Field<T> field, SmoothCurve<T> curve) {
		this.field = field;
		this.curve = curve;
		this.ring = new PolynomialRing<T>(this.field, this.curve.getEmbeddingDimension() + 1, Polynomial.LEX);
	}
	public FunctionField(Field<T> field, SmoothCurve<T> curve, PolynomialRing<T> ring) {
		this.field = field;
		this.curve = curve;
		this.ring = ring;
	}
	@Override
	public boolean isFinite() {
		return false;
	}
	@Override
	public RationalFunction<T> zero() {
		return new RationalFunction<T>(this.field, this.ring.zero(), this.ring.one(), this.curve, this.ring);
	}
	@Override
	public RationalFunction<T> one() {
		return new RationalFunction<T>(this.field, this.ring.one(), this.ring.one(), this.curve, this.ring);
	}
	@Override
	public int characteristic() {
		return this.field.characteristic();
	}
	@Override
	public RationalFunction<T> add(RationalFunction<T> t1,
			RationalFunction<T> t2) {
		Polynomial<T> numerator = this.ring.add(this.ring.multiply(t1.getNumerator(), t2.getDenominator()), this.ring.multiply(t2.getNumerator(), t1.getDenominator()));
		Polynomial<T> denominator = this.ring.multiply(t1.getDenominator(), t2.getDenominator());
		return new RationalFunction<T>(this.field, numerator, denominator, this.curve, this.ring);
	}
	@Override
	public RationalFunction<T> negative(RationalFunction<T> t) {
		return new RationalFunction<T>(this.field, this.ring.negative(t.getNumerator()), t.getDenominator(), this.curve, this.ring);
	}
	@Override
	public RationalFunction<T> multiply(RationalFunction<T> t1,
			RationalFunction<T> t2) {
		Polynomial<T> numerator = this.ring.multiply(t1.getNumerator(), t2.getNumerator());
		Polynomial<T> denominator = this.ring.multiply(t1.getDenominator(), t2.getDenominator());
		return new RationalFunction<T>(this.field, numerator, denominator, this.curve, this.ring);
	}
	@Override
	public RationalFunction<T> inverse(RationalFunction<T> t) {
		return new RationalFunction<T>(this.field, t.getDenominator(), t.getNumerator(), this.curve, this.ring);
	}
	@Override
	public RationalFunction<T> getRandomElement() {
		return null;/*
		Polynomial<T> nonzero = this.ring.getRandomElement();
		while (nonzero.equals(this.ring.zero())) {
			nonzero = this.ring.getRandomElement();
		}
		return new RationalFunction<T>(field.getField(), this.ring.getRandomElement(), nonzero);*/
	}
	@Override
	public int getNumberOfElements() {
		return -1;
	}
	@Override
	public Iterable<RationalFunction<T>> getElements() throws InfinityException {
		throw new InfinityException();
	}
}
