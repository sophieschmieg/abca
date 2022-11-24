package cryptography;

import java.io.IOException;
import java.io.InputStream;
import java.util.Arrays;

import cryptography.interfaces.ExtendedOutputFunction;
import cryptography.interfaces.HashFunction;
import cryptography.interfaces.PseudoRandomFunction;
import util.MiscAlgorithms;

public class Sha3 implements ExtendedOutputFunction, PseudoRandomFunction, HashFunction {
	private int width;
	private int capacity;
	private int length;
	private int l;
	private int w;

	public static final HashFunction SHA3_224 = new Sha3(6, 448, 224);
	public static final HashFunction SHA3_256 = new Sha3(6, 512, 256);
	public static final HashFunction SHA3_384 = new Sha3(6, 768, 384);
	public static final HashFunction SHA3_512 = new Sha3(6, 1024, 512);
	public static final ExtendedOutputFunction SHAKE_128 = new Sha3(6, 256, -1);
	public static final ExtendedOutputFunction SHAKE_256 = new Sha3(6, 512, -1);
	public static final PseudoRandomFunction SHAKE_128_PRF = new Sha3(6, 256, -1);
	public static final PseudoRandomFunction SHAKE_256_PRF = new Sha3(6, 512, -1);

	private Sha3(int l, int capacity, int length) {
		this.l = l;
		if (l < 3 || l > 6) {
			throw new ArithmeticException("wrong l");
		}
		this.w = 1 << l;
		this.width = 25 * w;
		this.capacity = capacity;
		this.length = length > 0 ? length / 8 : length;
	}

	@Override
	public InputStream stream(byte[] key, byte[] input) {
		byte[] concat = Arrays.copyOf(key, key.length + input.length);
		System.arraycopy(input, 0, concat, key.length, input.length);
		return stream(concat);
	}

	@Override
	public InputStream stream(byte[] input) {
		return keccak(input, 0x0f, 4);
	}

	@Override
	public byte[] evaluate(byte[] input) {
		InputStream in = keccak(input, 2, 2);
		byte[] result = new byte[length];
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

	private boolean bit(byte[] input, int index) {
		int byteIndex = index / 8;
		int bitIndex = index % 8;
		return (input[byteIndex] & (1 << bitIndex)) != 0;
	}

	private void setBit(byte[] data, boolean value, int index) {
		int byteIndex = index / 8;
		int bitIndex = index % 8;
		if (value) {
			data[byteIndex] |= (1 << bitIndex);
		} else {
			data[byteIndex] &= ~(1 << bitIndex);
		}
	}

	private boolean matrix(byte[] input, int x, int y, int z) {
		return bit(input, index(x, y, z));
	}

	private void setMatrix(byte[] input, boolean value, int x, int y, int z) {
		setBit(input, value, index(x, y, z));
	}

	private int index(int x, int y, int z) {
		return w * (5 * y + x) + z;
	}

	private void theta(byte[] state) {
		boolean[][] c = new boolean[5][w];
		for (int x = 0; x < 5; x++) {
			for (int z = 0; z < w; z++) {
				boolean xor = false;
				for (int y = 0; y < 5; y++) {
					xor ^= matrix(state, x, y, z);
				}
				c[x][z] = xor;
			}
		}
		for (int x = 0; x < 5; x++) {
			for (int z = 0; z < w; z++) {
				boolean d = c[(x + 4) % 5][z] ^ c[(x + 1) % 5][(z - 1 + w) % w];
				for (int y = 0; y < 5; y++) {
					setMatrix(state, matrix(state, x, y, z) ^ d, x, y, z);
				}
			}
		}
	}

	private void rho(byte[] state) {
		int x = 1;
		int y = 0;
		byte[] stateCopy = Arrays.copyOf(state, state.length);
		for (int t = 0; t < 24; t++) {
			int offset = ((t + 1) * (t + 2) / 2) % w;
			offset = (w - offset) % w;
			for (int z = 0; z < w; z++) {
				setMatrix(state, matrix(stateCopy, x, y, (z + offset) % w), x, y, z);
			}
			int tmp = x;
			x = y;
			y = (2 * tmp + 3 * y) % 5;
		}
	}

	private void pi(byte[] state) {
		byte[] stateCopy = Arrays.copyOf(state, state.length);
		for (int x = 0; x < 5; x++) {
			for (int y = 0; y < 5; y++) {
				for (int z = 0; z < w; z++) {
					setMatrix(state, matrix(stateCopy, (x + 3 * y) % 5, x, z), x, y, z);
				}
			}
		}
	}

	private void chi(byte[] state) {
		byte[] stateCopy = Arrays.copyOf(state, state.length);
		for (int x = 0; x < 5; x++) {
			for (int y = 0; y < 5; y++) {
				for (int z = 0; z < w; z++) {
					setMatrix(state, matrix(stateCopy, x, y, z)
							^ ((matrix(stateCopy, (x + 1) % 5, y, z) ^ true) && matrix(stateCopy, (x + 2) % 5, y, z)),
							x, y, z);
				}
			}
		}
	}

	private boolean rc(int t) {
		int r = 0x01;
		t %= 255;
		for (int i = 0; i < t; i++) {
			r <<= 1;
			int mask = -(r >> 8) & 0x71;
			r ^= mask;
			r &= 0xff;
		}
		return (r & 1) != 0;
	}

	private void iota(byte[] state, int roundIndex) {
		for (int j = 0; j <= l; j++) {
			int z = (1 << j) - 1;
			setMatrix(state, matrix(state, 0, 0, z) ^ rc(j + 7 * roundIndex), 0, 0, z);
		}
	}

	private void round(byte[] state, int roundIndex) {
		theta(state);
		rho(state);
		pi(state);
		chi(state);
		iota(state, roundIndex);
	}

	private void keccakP(byte[] state, int numberOfRounds) {
		for (int roundIndex = 12 + 2 * l - numberOfRounds; roundIndex < 12 + 2 * l; roundIndex++) {
			round(state, roundIndex);
		}
	}

	private void keccakF(byte[] state) {
		keccakP(state, 12 + 2 * l);
	}

	private InputStream keccak(byte[] input, int additionalPadding, int additionalPaddingLength) {
		return new KeccakInputStream(input, additionalPadding, additionalPaddingLength);
	}

	private class KeccakInputStream extends InputStream {
		private byte[] state;
		private int index;
		private int rate;

		private KeccakInputStream(byte[] input, int additionalPadding, int additionalPaddingLength) {
			this.rate = width - capacity;
			if (rate % 8 != 0) {
				throw new ArithmeticException("Capacity wrong!");
			}
			int paddingLength = (8 * input.length % rate);
			paddingLength = (-paddingLength + rate - 2) % rate;
			byte[] padded = Arrays.copyOf(input, input.length + MiscAlgorithms.DivRoundUp(paddingLength + 2, 8));
			padded[input.length] |= additionalPadding;
			padded[input.length] |= 1 << additionalPaddingLength;
			padded[padded.length - 1] |= 0x80;
			int numAbsorbtions = 8 * padded.length / rate;
			this.state = new byte[width / 8];
			for (int i = 0; i < numAbsorbtions; i++) {
				for (int j = 0; j < rate / 8; j++) {
					state[j] ^= padded[i * rate / 8 + j];
				}
				keccakF(state);
			}
			index = 0;
		}

		@Override
		public int read() throws IOException {
			if (index >= rate / 8) {
				keccakF(state);
				index = 0;
			}
			index++;
			return ((int) state[index - 1]) & 0xff;
		}
	}

	@Override
	public int outputLength() {
		return length;
	}

	@Override
	public int blockSize() {
		return capacity / 8;
	}
}
