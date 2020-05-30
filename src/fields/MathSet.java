package fields;

import java.math.BigInteger;

public interface MathSet<T extends Element> extends Iterable<T> {
	public T getRandomElement();
	public boolean isFinite();
	public BigInteger getNumberOfElements() throws InfinityException;
}
