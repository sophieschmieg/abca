package varieties;

import java.math.BigInteger;
import java.util.Iterator;

import fields.exceptions.InfinityException;
import fields.helper.AbstractField;
import fields.helper.CoordinateRing;
import fields.helper.CoordinateRing.CoordinateRingElement;
import fields.interfaces.Element;
import fields.interfaces.Field;
import fields.interfaces.PolynomialRing;
import fields.polynomials.AbstractPolynomialRing;
import fields.polynomials.Monomial;
import varieties.curves.SmoothCurve;

public class FunctionField<T extends Element<T>> extends AbstractField<RationalFunction<T>> {
	private Field<T> field;
	private PolynomialRing<T> ring;
	private CoordinateRing<T> coordRing;
	private SmoothCurve<T> curve;

	public FunctionField(Field<T> field, SmoothCurve<T> curve) {
		this(field, curve, AbstractPolynomialRing.getPolynomialRing(field, curve.getEmbeddingDimension() + 1, Monomial.GREVLEX));
	}

	public FunctionField(Field<T> field, SmoothCurve<T> curve, PolynomialRing<T> ring) {
		this.field = field;
		this.curve = curve;
		this.ring = ring;
		this.coordRing = this.curve.getCoordinateRing();
	}

	@Override
	public boolean isFinite() {
		return false;
	}

	public RationalFunction<T> getEmbedding(CoordinateRingElement<T> t) {
		return getFunction(t, coordRing.one());
	}

	public RationalFunction<T> getEmbedding(T t) {
		return getEmbedding(coordRing.getEmbedding(t));
	}

	public RationalFunction<T> getFunction(CoordinateRingElement<T> numerator, CoordinateRingElement<T> denominator) {
		return new RationalFunction<T>(field, numerator, denominator, curve, ring);
	}

	@Override
	public RationalFunction<T> zero() {
		return getEmbedding(field.zero());
	}

	@Override
	public RationalFunction<T> one() {
		return getEmbedding(field.one());
	}

	@Override
	public BigInteger characteristic() {
		return this.field.characteristic();
	}

	@Override
	public RationalFunction<T> add(RationalFunction<T> t1, RationalFunction<T> t2) {
		CoordinateRingElement<T> numerator = coordRing.add(coordRing.multiply(t1.getNumerator(), t2.getDenominator()),
				coordRing.multiply(t2.getNumerator(), t1.getDenominator()));
		CoordinateRingElement<T> denominator = coordRing.multiply(t1.getDenominator(), t2.getDenominator());
		RationalFunction<T> rf = new RationalFunction<T>(this.field, numerator, denominator, this.curve, this.ring);
		rf.simplify();
		return rf;
	}

	@Override
	public RationalFunction<T> negative(RationalFunction<T> t) {
		return new RationalFunction<T>(this.field, coordRing.negative(t.getNumerator()), t.getDenominator(), this.curve,
				this.ring);
	}

	@Override
	public RationalFunction<T> multiply(RationalFunction<T> t1, RationalFunction<T> t2) {
		CoordinateRingElement<T> numerator = coordRing.multiply(t1.getNumerator(), t2.getNumerator());
		CoordinateRingElement<T> denominator = coordRing.multiply(t1.getDenominator(), t2.getDenominator());
		RationalFunction<T> rf = new RationalFunction<T>(this.field, numerator, denominator, this.curve, this.ring);
		rf.simplify();
		return rf;
	}

	@Override
	public RationalFunction<T> inverse(RationalFunction<T> t) {
		return new RationalFunction<T>(this.field, t.getDenominator(), t.getNumerator(), this.curve, this.ring);
	}

	@Override
	public RationalFunction<T> getRandomElement() {
		return null;
	}

	@Override
	public BigInteger getNumberOfElements() {
		return BigInteger.valueOf(-1);
	}

	@Override
	public Iterator<RationalFunction<T>> iterator() {
		throw new InfinityException();
	}
}
