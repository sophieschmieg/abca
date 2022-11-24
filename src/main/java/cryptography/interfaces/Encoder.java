package cryptography.interfaces;

import java.util.function.Function;

import fields.interfaces.Element;

public interface Encoder<T extends Element<T>> {
	public static <T extends Element<T>> Encoder<T> fromFunction(Function<T, byte[]> function) {
		return new Encoder<>() {
			@Override
			public byte[] encode(T t) {
				return function.apply(t);
			}
		};
	}

	byte[] encode(T t);
}
