package fields.interfaces;

import java.io.IOException;
import java.io.StringReader;
import java.math.BigInteger;

import fields.exceptions.InfinityException;
import util.PeekableReader;

public interface MathSet<T extends Element<T>> extends Iterable<T> {
	public enum Exactness {
		EXACT, FLOATING_POINT, FIXED_POINT
	}

	public Exactness exactness();

	public T getRandomElement();

	public default T getRandomElement(T max) {
		return getRandomElement();
	}

	public boolean isFinite();

	public BigInteger getNumberOfElements() throws InfinityException;

	default public T parse(PeekableReader reader) throws IOException {
		throw new IOException("Not implemented for " + this);
	}

	default public T parse(String text) throws IOException {
		PeekableReader reader = new PeekableReader(new StringReader(text));
		T result = parse(reader);
		if (reader.peek() >= 0) {
			throw new IOException("Parser did not consume full element!");
		}
		return result;
	}
}
