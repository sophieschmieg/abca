package cryptography;

import java.security.MessageDigest;
import java.security.NoSuchAlgorithmException;

import cryptography.interfaces.HashFunction;

public class Sha2 implements HashFunction {
	public static final HashFunction SHA2_256 = new Sha2("SHA-256", 32, 64);
	public static final HashFunction SHA2_512 = new Sha2("SHA-512", 64, 128);
	
	private String algorithm;
	private int outputLength;
	private int blockSize;

	private Sha2(String algorithm, int outputLength, int blockSize) {
		this.algorithm = algorithm;
		this.outputLength = outputLength;
		this.blockSize = blockSize;
	}

	@Override
	public byte[] evaluate(byte[] input) {
		try {
			MessageDigest digest = MessageDigest.getInstance(algorithm);
			return digest.digest(input);
		} catch (NoSuchAlgorithmException e) {
			throw new RuntimeException("Algorithm not found", e);
		}
	}

	@Override
	public int outputLength() {
		return outputLength;
	}
	
	@Override
	public int blockSize() {
		return blockSize;
	}
}
