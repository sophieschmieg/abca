package cryptography;

import java.io.ByteArrayInputStream;
import java.io.InputStream;

import cryptography.interfaces.HashFunction;
import cryptography.interfaces.PseudoRandomFunction;

public class Hmac implements PseudoRandomFunction {
	private HashFunction hash;
	private int blockSize;
	private final static byte outerPad = 0x5c;
	private final static byte innerPad = 0x36;

	public Hmac(HashFunction hash) {
		this.hash = hash;
		this.blockSize = hash.blockSize();
	}

	@Override
	public InputStream stream(byte[] key, byte[] input) {
		byte[] outerKey = new byte[blockSize];
		byte[] innerKey = new byte[blockSize];
		if (key.length > blockSize) {
			key = hash.evaluate(key);
		}
		for (int i = 0; i < blockSize; i++) {
			outerKey[i] = (byte) ((i < key.length ? key[i] : (byte) 0x00) ^ outerPad);
			innerKey[i] = (byte) ((i < key.length ? key[i] : (byte) 0x00) ^ innerPad);
		}
		byte[] innerInput = new byte[blockSize + input.length];
		System.arraycopy(innerKey, 0, innerInput, 0, blockSize);
		System.arraycopy(input, 0, innerInput, blockSize, input.length);
		byte[] outerData = hash.evaluate(innerInput);
		byte[] outerInput = new byte[blockSize + outerData.length];
		System.arraycopy(outerKey, 0, outerInput, 0, blockSize);
		System.arraycopy(outerData, 0, outerInput, blockSize, outerData.length);
		byte[] result = hash.evaluate(outerInput);
		return new ByteArrayInputStream(result);
	}

	@Override
	public int outputLength() {
		return hash.outputLength();
	}
}
