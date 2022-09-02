package cryptography;

import java.io.IOException;
import java.io.InputStream;
import java.util.Arrays;

import cryptography.interfaces.ExtendedOutputFunction;
import fields.helper.AbstractElement;

public class VariableLengthKey extends AbstractElement<VariableLengthKey> {
	private InputStream in;
	private byte[] seed;

	public VariableLengthKey(byte[] seed, ExtendedOutputFunction xof) {
		this.in = xof.stream(seed);
		this.seed = Arrays.copyOf(seed, seed.length);
	}

	public byte[] key(int length) {
		byte[] key = new byte[length];
		for (int i = 0; i < length; i++) {
			try {
				int read = in.read();
				if (read < 0) {
					throw new RuntimeException("Inputstream reached end of stream!");
				}
				key[i] = (byte) read;
			} catch (IOException e) {
				throw new RuntimeException("Inputstream failed", e);
			}
		}
		return key;
	}

	@Override
	public int compareTo(VariableLengthKey o) {
		return Arrays.compare(seed, o.seed);
	}
}
