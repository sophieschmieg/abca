package fields.interfaces;

import java.math.BigInteger;

import fields.exceptions.InfinityException;

public interface MathSet<T extends Element<T>> extends Iterable<T> {
	public T getRandomElement();
	public boolean isFinite();
	public BigInteger getNumberOfElements() throws InfinityException;
}
