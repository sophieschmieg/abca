package cryptography;

import java.util.Arrays;

import cryptography.interfaces.HashFunction;

public class TruncatedHashFunction implements HashFunction {
	private HashFunction hash;
	private int length;

	public TruncatedHashFunction(HashFunction hash, int length) {
		this.hash = hash;
		this.length = length;
	}

	@Override
	public byte[] evaluate(byte[] input) {
		return Arrays.copyOf(hash.evaluate(input), length);
	}
	
	@Override
	public int blockSize() {
		return hash.blockSize();
	}
	
	@Override
	public int outputLength() {
		return length;
	}
}
