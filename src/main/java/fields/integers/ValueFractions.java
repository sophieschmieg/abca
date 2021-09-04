package fields.integers;

import java.math.BigInteger;
import java.util.Iterator;

import fields.exceptions.InfinityException;
import fields.floatingpoint.Reals;
import fields.floatingpoint.Reals.Real;
import fields.helper.AbstractField;
import fields.integers.Rationals.Fraction;
import fields.interfaces.ValueField;

public class ValueFractions extends AbstractField<Fraction> implements ValueField<Fraction> {
	private Reals r = Reals.r(1024);
	
	private Rationals q = Rationals.q();

	ValueFractions() {
	}

	@Override
	public boolean isComplete() {
		return false;
	}

	@Override
	public Real value(Fraction t) {
		return r.abs(r.getFraction(t));
		//double absNumerator = t.getNumerator().getValue().abs().doubleValue();
		//double denominator = t.getDenominator().getValue().doubleValue();
		//return absNumerator / denominator;
	}
	
	public Fraction abs(Fraction t) {
		if (t.compareTo(zero()) >= 0) {
			return t;
		}
		return negative(t);
	}

	@Override
	public Fraction zero() {
		return q.zero();
	}

	@Override
	public Fraction one() {
		return q.one();
	}

	@Override
	public BigInteger characteristic() {
		return BigInteger.ZERO;
	}

	@Override
	public Fraction add(Fraction t1, Fraction t2) {
		return q.add(t1, t2);
	}

	@Override
	public Fraction negative(Fraction t) {
		return q.negative(t);
	}

	@Override
	public Fraction multiply(Fraction t1, Fraction t2) {
		return q.multiply(t1, t2);
	}

	@Override
	public Fraction inverse(Fraction t) {
		return q.inverse(t);
	}

	@Override
	public Exactness exactness() {
		return Exactness.EXACT;
	}

	@Override
	public Fraction getRandomElement() {
		return q.getRandomElement();
	}

	@Override
	public boolean isFinite() {
		return false;
	}

	@Override
	public BigInteger getNumberOfElements() throws InfinityException {
		throw new InfinityException();
	}

	@Override
	public Iterator<Fraction> iterator() {
		return q.iterator();
	}

	public Fraction round(Fraction t) {
		return q.getInteger(t.round());
	}

}
