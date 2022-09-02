package cryptography.interfaces;

import java.io.IOException;
import java.io.InputStream;

public interface ExtendedOutputFunction {
	InputStream stream(byte[] input);

	default byte[] evaluate(byte[] input, int length) {
		byte[] result = new byte[length];
		InputStream in = stream(input);
		for (int i = 0; i < length; i++) {
			try {
				int read = in.read();
				if (read < 0) {
					throw new ArithmeticException("XOF function reached end of stream!");
				}
				result[i] = (byte) read;
			} catch (IOException e) {
				throw new RuntimeException("XOF function had an IOException", e);
			}
		}
		return result;
	}
}
