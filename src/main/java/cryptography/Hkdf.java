package cryptography;

import java.io.IOException;
import java.io.InputStream;

import cryptography.interfaces.ExtendedOutputFunction;
import cryptography.interfaces.HashFunction;
import cryptography.interfaces.PseudoRandomFunction;

public class Hkdf implements ExtendedOutputFunction {
	private PseudoRandomFunction prf;
	private byte[] salt;

	public Hkdf(HashFunction hash) {
		this(new Hmac(hash));
	}

	public Hkdf(HashFunction hash, byte[] salt) {
		this(new Hmac(hash), salt);
	}

	public Hkdf(PseudoRandomFunction prf) {
		this(prf, new byte[prf.outputLength()]);
	}

	public Hkdf(PseudoRandomFunction prf, byte[] salt) {
		this.prf = prf;
		this.salt = salt;
	}

	private byte[] extract(byte[] input) {
		return prf.evaluate(salt, input, prf.outputLength());
	}

	private InputStream expand(byte[] ikm) {
		return new InputStream() {
			private int counter = 0;
			private int index = 0;
			private byte[] currentBlock = new byte[0];

			@Override
			public int read() throws IOException {
				if (index >= currentBlock.length) {
					if (counter > 255) {
						return -1;
					}
					byte[] input = new byte[currentBlock.length + 1];
					System.arraycopy(currentBlock, 0, input, 0, currentBlock.length);
					input[currentBlock.length] = (byte) counter;
					counter++;
					index = 0;
					currentBlock = prf.evaluate(ikm, input, prf.outputLength());
				}
				index++;
				return currentBlock[index - 1] & 0xff;
			}

		};
	}

	@Override
	public InputStream stream(byte[] input) {
		return expand(extract(input));
	}

}
