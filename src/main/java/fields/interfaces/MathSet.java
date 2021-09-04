package fields.interfaces;

import java.math.BigInteger;

import fields.exceptions.InfinityException;

public interface MathSet<T extends Element<T>> extends Iterable<T> {
	public enum Exactness{
		EXACT,
		FLOATING_POINT,
		FIXED_POINT
	}
	public Exactness exactness();
	public T getRandomElement();
	public boolean isFinite();
	public BigInteger getNumberOfElements() throws InfinityException;
}
