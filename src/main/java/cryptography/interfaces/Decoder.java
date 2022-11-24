package cryptography.interfaces;

import java.util.function.Function;

import fields.interfaces.Element;

public interface Decoder<T extends Element<T>> {
	public static <T extends Element<T>> Decoder<T> fromFunction(Function<byte[], T> function) {
		return new Decoder<>() {
			@Override
			public T decode(byte[] t) {
				return function.apply(t);
			}
		};
	}

	T decode(byte[] t);
}
