package cryptography;

import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.security.SecureRandom;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;
import java.util.Optional;

import cryptography.Dilithium.DilithiumPrivateKey;
import cryptography.Dilithium.DilithiumPublicKey;
import cryptography.Dilithium.DilithiumSignature;
import cryptography.interfaces.ExtendedOutputFunction;
import cryptography.interfaces.SignatureScheme;
import fields.finitefields.PrimeField;
import fields.finitefields.PrimeField.PFE;
import fields.helper.AbstractElement;
import fields.helper.NTTRing;
import fields.helper.NTTRing.NTT;
import fields.integers.Integers;
import fields.integers.Integers.IntE;
import fields.interfaces.UnivariatePolynomial;
import fields.vectors.FreeModule;
import fields.vectors.Matrix;
import fields.vectors.MatrixModule;
import fields.vectors.Vector;
import util.MiscAlgorithms;
import util.Pair;

public class Dilithium
		implements SignatureScheme<DilithiumSignature, DilithiumPublicKey, DilithiumPrivateKey, Dilithium> {
	public static class DilithiumPublicKey extends AbstractElement<DilithiumPublicKey> {
		private byte[] rho;
		private Vector<IntE> t1;

		private DilithiumPublicKey(byte[] rho, Vector<IntE> t1) {
			super();
			this.rho = rho;
			this.t1 = t1;
		}

		public byte[] getRho() {
			return rho;
		}

		public Vector<IntE> getT1() {
			return t1;
		}

		@Override
		public int compareTo(DilithiumPublicKey o) {
			int cmp = Arrays.compare(rho, o.rho);
			if (cmp != 0) {
				return cmp;
			}
			return t1.compareTo(o.t1);
		}

	}

	public static class DilithiumPrivateKey extends AbstractElement<DilithiumPrivateKey> {
		private DilithiumPublicKey publicKey;
		private byte[] rho;
		private byte[] k;
		private byte[] tr;
		private Vector<NTT<PFE>> s;
		private Vector<NTT<PFE>> e;
		private Vector<NTT<PFE>> t0;

		public DilithiumPrivateKey(DilithiumPublicKey publicKey, byte[] k, byte[] tr, Vector<NTT<PFE>> s,
				Vector<NTT<PFE>> e, Vector<NTT<PFE>> t0) {
			this.publicKey = publicKey;
			this.rho = publicKey.getRho();
			this.k = k;
			this.tr = tr;
			this.s = s;
			this.e = e;
			this.t0 = t0;
		}

		public DilithiumPublicKey getPublicKey() {
			return publicKey;
		}

		public byte[] getRho() {
			return rho;
		}

		public byte[] getK() {
			return k;
		}

		public byte[] getTr() {
			return tr;
		}

		public Vector<NTT<PFE>> getS() {
			return s;
		}

		public Vector<NTT<PFE>> getE() {
			return e;
		}

		public Vector<NTT<PFE>> getT0() {
			return t0;
		}

		@Override
		public int compareTo(DilithiumPrivateKey o) {
			return publicKey.compareTo(o.publicKey);
		}

	}

	public static class DilithiumSignature extends AbstractElement<DilithiumSignature> {
		private Vector<NTT<PFE>> z;
		private List<Boolean> hints;
		private byte[] cSeed;

		private DilithiumSignature(Vector<NTT<PFE>> z, List<Boolean> hints, byte[] cSeed) {
			this.z = z;
			this.hints = hints;
			this.cSeed = cSeed;
		}

		public Vector<NTT<PFE>> getZ() {
			return z;
		}

		public List<Boolean> getHints() {
			return hints;
		}

		public byte[] getCSeed() {
			return cSeed;
		}

		@Override
		public int compareTo(DilithiumSignature o) {
			int cmp = z.compareTo(o.z);
			if (cmp != 0) {
				return cmp;
			}
			for (int i = 0; i < hints.size(); i++) {
				if (hints.get(i) && !o.hints.get(i)) {
					return 1;
				}
				if (!hints.get(i) && hints.get(i)) {
					return -1;
				}
			}
			return Arrays.compare(cSeed, o.cSeed);
		}

	}

	private IntE q;
	private int logQ;
	private PrimeField field;
	private int degree;
	private int rows;
	private int cols;
	private int d;
	private int tau;
	private int eta;
	private int gamma1;
	private int log2Gamma1;
	private int gamma2;
	private int omega;
	private int w1bits;
	private int lambda;
	private NTTRing<PFE, PrimeField> cylcotomic;
	private FreeModule<NTT<PFE>> freeModule;
	private MatrixModule<NTT<PFE>> matrixModule;
	private ExtendedOutputFunction hashG;
	private ExtendedOutputFunction hashH;
	// private ExtendedOutputFunction crh;

	public final static int PRIME = (1 << 23) - (1 << 13) + 1;

	public final static Dilithium DILITHIUM2 = new Dilithium(4, 4, 2, 13, 39, 17, (PRIME - 1) / 88, 80, 128);
	public final static Dilithium DILITHIUM3 = new Dilithium(6, 5, 4, 13, 49, 19, (PRIME - 1) / 32, 55, 192);
	public final static Dilithium DILITHIUM5 = new Dilithium(8, 7, 2, 13, 60, 19, (PRIME - 1) / 32, 75, 256);

	public Dilithium(int k, int l, int eta, int d, int tau, int log2Gamma1, int gamma2, int omega, int lambda) {
		this(PRIME, 256, 1753, k, l, eta, d, tau, log2Gamma1, gamma2, omega, lambda, Sha3.SHAKE_256, Sha3.SHAKE_128);
	}

	public Dilithium(int q, int degree, int r, int k, int l, int eta, int d, int tau, int log2Gamma1, int gamma2,
			int omega, int lambda, ExtendedOutputFunction hashH, ExtendedOutputFunction hashG) {
		Integers z = Integers.z();
		this.q = z.getInteger(q);
		this.degree = degree;
		this.rows = k;
		this.cols = l;
		this.eta = eta;
		this.field = PrimeField.getPrimeField(q);
		logQ = 0;
		int power = 1;
		while (power < q) {
			power <<= 1;
			logQ++;
		}
		this.cylcotomic = new NTTRing<PrimeField.PFE, PrimeField>(field, degree, field.getInteger(r));
		this.matrixModule = new MatrixModule<>(cylcotomic, rows, cols);
		this.freeModule = matrixModule.codomain();
		this.d = d;
		this.tau = tau;
		this.log2Gamma1 = log2Gamma1;
		this.gamma1 = 1 << log2Gamma1;
		this.gamma2 = gamma2;
		this.omega = omega;
		this.w1bits = MiscAlgorithms.bitLength((q - 1) / (2 * gamma2) - 1);
		this.lambda = lambda;
		this.hashH = hashH;
		this.hashG = hashG;
	}

	public int getDegree() {
		return degree;
	}

	public int getRows() {
		return rows;
	}

	public int getCols() {
		return cols;
	}

	public int getEta() {
		return eta;
	}

	public NTTRing<PFE, PrimeField> getNTTRing() {
		return cylcotomic;
	}

	public MatrixModule<NTT<PFE>> getMatrixModule() {
		return matrixModule;
	}

	public FreeModule<NTT<PFE>> getFreeModule() {
		return freeModule;
	}

	private Pair<Vector<IntE>, Vector<NTT<PFE>>> power2Round(Vector<NTT<PFE>> t, int exp) {
		List<IntE> result1 = new ArrayList<>();
		List<NTT<PFE>> result2 = new ArrayList<>();
		for (int i = 0; i < t.dimension(); i++) {
			UnivariatePolynomial<PFE> coeff = cylcotomic.asPolynomial(t.get(i + 1));
			List<PFE> coeff2 = new ArrayList<>();
			for (int j = 0; j < degree; j++) {
				Pair<IntE, IntE> power2Round = power2Round(field.liftToInteger(coeff.univariateCoefficient(j)), exp);
				result1.add(power2Round.getFirst());
				coeff2.add(field.reduce(power2Round.getSecond()));
			}
			result2.add(cylcotomic.fromPolynomial(field.getUnivariatePolynomialRing().getPolynomial(coeff2)));
		}
		return new Pair<>(new Vector<>(result1), new Vector<>(result2));
	}

	private Pair<IntE, IntE> power2Round(IntE t, int exp) {
		Integers z = Integers.z();
		t = z.remainder(t, q);
		IntE power = z.power(z.getInteger(2), exp);
		IntE halfPower = z.divideChecked(power, z.getInteger(2));
		IntE r0 = z.remainder(t, power);
		if (r0.compareTo(halfPower) > 0) {
			r0 = z.subtract(r0, power);
		}
		IntE quotient = z.divideChecked(z.subtract(t, r0), power);
		return new Pair<>(quotient, r0);
	}

	private Vector<NTT<PFE>> power2Inflate(Vector<IntE> t, int exp) {
		Integers z = Integers.z();
		List<NTT<PFE>> result = new ArrayList<>();
		for (int i = 0; i < t.dimension() / degree; i++) {
			List<PFE> coeff = new ArrayList<>();
			for (int j = 0; j < degree; j++) {
				coeff.add(field.getInteger(z.multiply(t.get((degree * i) + j + 1), z.power(z.getInteger(2), exp))));
			}
			result.add(cylcotomic.fromPolynomial(field.getUnivariatePolynomialRing().getPolynomial(coeff)));
		}
		return new Vector<>(result);
	}

	private Pair<Integer, Integer> decompose(int t, int alpha) {
		t %= q.intValueExact();
		int remainder = t % alpha;
		if (2 * remainder > alpha) {
			remainder -= alpha;
		}
		if (t - remainder == q.intValueExact() - 1) {
			return new Pair<>(0, remainder - 1);
		}
		return new Pair<>((t - remainder) / alpha, remainder);
	}

	private Pair<Vector<IntE>, Vector<NTT<PFE>>> decompose(Vector<NTT<PFE>> t, int alpha) {
		Integers z = Integers.z();
		List<IntE> result1 = new ArrayList<>();
		List<NTT<PFE>> result2 = new ArrayList<>();
		for (int i = 0; i < t.dimension(); i++) {
			UnivariatePolynomial<PFE> coeff = cylcotomic.asPolynomial(t.get(i + 1));
			List<PFE> resultCoeff = new ArrayList<>();
			for (int j = 0; j < degree; j++) {
				Pair<Integer, Integer> decomposed = decompose(coeff.univariateCoefficient(j).getValue().intValueExact(),
						alpha);
				result1.add(z.getInteger(decomposed.getFirst()));
				resultCoeff.add(field.getInteger(decomposed.getSecond()));
			}
			result2.add(cylcotomic.fromPolynomial(field.getUnivariatePolynomialRing().getPolynomial(resultCoeff)));
		}
		return new Pair<>(new Vector<>(result1), new Vector<>(result2));
	}

	private int highBits(int t, int alpha) {
		return decompose(t, alpha).getFirst();
	}

	private Vector<IntE> highBits(Vector<NTT<PFE>> t, int alpha) {
		return decompose(t, alpha).getFirst();
	}

	private Vector<NTT<PFE>> lowBits(Vector<NTT<PFE>> t, int alpha) {
		return decompose(t, alpha).getSecond();
	}

	private static class IntPacker {
		private ByteArrayOutputStream result;
		private int current;
		private int currentOffset;
		private int log;
		private boolean closed;

		public IntPacker(int log) {
			this.log = log;
			this.current = 0;
			this.currentOffset = 0;
			this.result = new ByteArrayOutputStream();
			this.closed = false;
		}

		public void write(int t) {
			if (closed) {
				throw new RuntimeException("IntPacker closed!");
			}
			if (t < 0) {
				throw new RuntimeException("can only write non-negative integers");
			}
			if (t >= (1 << log)) {
				throw new RuntimeException("can only write integers with " + log + " bits");
			}
			int logToWrite = log;
			while (logToWrite > 0) {
				int bits = Math.min(logToWrite, 8 - currentOffset);
				int mod = t % (1 << bits);
				current |= mod << currentOffset;
				currentOffset += bits;
				logToWrite -= bits;
				t >>= bits;
				if (currentOffset == 8) {
					result.write(current);
					currentOffset = 0;
					current = 0;
				}
			}
		}

		public byte[] close() {
			if (currentOffset != 0) {
				throw new RuntimeException("Dangling bits in IntPacker!");
			}
			closed = true;
			try {
				result.close();
			} catch (IOException e) {
				throw new RuntimeException("ByteArrayOutputStream could not close!", e);
			}
			return result.toByteArray();

		}
	}

	private byte[] packVector(Vector<IntE> t, int bits) {
		IntPacker result = new IntPacker(bits);
		for (int i = 0; i < t.dimension(); i++) {
			result.write(t.get(i + 1).intValueExact());
		}
		return result.close();
	}

	private int infinityNorm(Vector<NTT<PFE>> t) {
		int result = 0;
		for (int i = 0; i < t.dimension(); i++) {
			UnivariatePolynomial<PFE> coeff = cylcotomic.asPolynomial(t.get(i + 1));
			for (int j = 0; j < degree; j++) {
				int value = coeff.univariateCoefficient(j).getValue().intValueExact();
				if (PRIME - value < value) {
					value = PRIME - value;
				}
				if (value > result) {
					result = value;
				}
			}
		}
		return result;
	}

	private int oneNorm(List<Boolean> t) {
		int result = 0;
		for (boolean c : t) {
			if (c) {
				result++;
			}
		}
		return result;
	}

	private boolean makeHint(int z, int r, int alpha) {
		int highR = highBits(r, alpha);
		int highV = highBits(z + r, alpha);
		return highR != highV;
	}

	private List<Boolean> makeHints(Vector<NTT<PFE>> z, Vector<NTT<PFE>> r, int alpha) {
		List<Boolean> result = new ArrayList<>();
		for (int i = 0; i < z.dimension(); i++) {
			UnivariatePolynomial<PFE> zCoeff = cylcotomic.asPolynomial(z.get(i + 1));
			UnivariatePolynomial<PFE> rCoeff = cylcotomic.asPolynomial(r.get(i + 1));
			for (int j = 0; j < degree; j++) {
				result.add(makeHint(zCoeff.univariateCoefficient(j).getValue().intValueExact(),
						rCoeff.univariateCoefficient(j).getValue().intValueExact(), alpha));
			}
		}
		return result;
	}

	private int useHint(boolean hint, int r, int alpha) {
		int m = (q.intValueExact() - 1) / alpha;
		Pair<Integer, Integer> decomposed = decompose(r, alpha);
		if (hint && decomposed.getSecond() > 0) {
			return (decomposed.getFirst() + 1) % m;
		}
		if (hint && decomposed.getSecond() <= 0) {
			return (decomposed.getFirst() + m - 1) % m;
		}
		return decomposed.getFirst();
	}

	private Vector<IntE> useHints(List<Boolean> hints, Vector<NTT<PFE>> r, int alpha) {
		Integers z = Integers.z();
		List<IntE> result = new ArrayList<>();
		int index = 0;
		for (int i = 0; i < r.dimension(); i++) {
			UnivariatePolynomial<PFE> coeff = cylcotomic.asPolynomial(r.get(i + 1));
			for (int j = 0; j < degree; j++) {
				result.add(z.getInteger(
						useHint(hints.get(index), coeff.univariateCoefficient(j).getValue().intValueExact(), alpha)));
				index++;
			}
		}
		return new Vector<>(result);
	}

	private Iterator<IntE> reinterpretStreamUnbounded(InputStream stream) {
		return new Iterator<>() {
			private IntE next = null;
			private Integers z = Integers.z();

			@Override
			public boolean hasNext() {
				if (next != null) {
					return true;
				}
				try {
					int result;
					do {
						int b0 = stream.read();
						int b1 = stream.read();
						int b2 = stream.read();
						if (b0 < 0 || b1 < 0 || b2 < 0) {
							return false;
						}
						if (b2 > 127) {
							b2 -= 128;
						}
						result = b0 + 256 * b1 + 65536 * b2;
					} while (result >= q.intValueExact());
					next = z.getInteger(result);
				} catch (IOException e) {
					return false;
				}
				return true;
			}

			@Override
			public IntE next() {
				if (!hasNext()) {
					throw new ArithmeticException("Stream ended unexpectedly!");
				}
				IntE result = next;
				next = null;
				return result;
			}

		};
	}

	private Iterator<IntE> reinterpretStreamBounded(InputStream stream) {
		return new Iterator<>() {
			private IntE[] next = new IntE[2];
			private Integers z = Integers.z();

			private Optional<Integer> sample(int input) {
				if (eta == 2 && input < 15) {
					return Optional.of(2 - (input % 5));
				}
				if (eta == 4 && input < 9) {
					return Optional.of(4 - input);
				}
				return Optional.empty();
			}

			@Override
			public boolean hasNext() {
				if (next[0] != null) {
					return true;
				}
				try {
					do {
						int b = stream.read();
						if (b < 0) {
							return false;
						}
						int index = 0;
						Optional<Integer> low = sample(b % 16);
						if (low.isPresent()) {
							next[index] = z.getInteger(low.get());
							index++;
						}
						Optional<Integer> high = sample(b / 16);
						if (high.isPresent()) {
							next[index] = z.getInteger(high.get());
							index++;
						}

					} while (next[0] == null);
				} catch (IOException e) {
					return false;
				}
				return true;
			}

			@Override
			public IntE next() {
				if (!hasNext()) {
					throw new ArithmeticException("Stream ended unexpectedly!");
				}
				IntE result = next[0];
				next[0] = next[1];
				next[1] = null;
				return result;
			}

		};
	}
	/*
	 * private Iterator<IntE> reinterpretStream(int log, InputStream stream) {
	 * return new Iterator<>() { private Integers z = Integers.z(); private int bits
	 * = 0; private int leftOver = 0; private IntE next = null;
	 * 
	 * private boolean getBit() throws IOException { if (bits == 0) { leftOver =
	 * stream.read(); if (leftOver < 0) { throw new
	 * IOException("End of stream reached"); } bits = 8; } boolean result =
	 * (leftOver & 0x01) != 0; bits--; leftOver >>= 1; return result; }
	 * 
	 * @Override public boolean hasNext() { if (next != null) { return true; } next
	 * = z.zero(); IntE power = z.one(); try { for (int i = 0; i < log; i++) { if
	 * (getBit()) { next = z.add(next, power); } power = z.multiply(2, power); } }
	 * catch (IOException e) { return false; } return true; }
	 * 
	 * @Override public IntE next() { if (!hasNext()) { throw new
	 * ArithmeticException("Stream ended unexpectedly!"); } IntE result = next; next
	 * = null; return result; } }; }
	 */

	private NTT<PFE> sampleBounded(InputStream stream) {
		List<PFE> coeffs = new ArrayList<>();
		Iterator<IntE> it = reinterpretStreamBounded(stream);
		while (coeffs.size() < degree) {
			coeffs.add(field.getInteger(it.next()));
		}
		return cylcotomic.fromPolynomial(field.getUnivariatePolynomialRing().getPolynomial(coeffs));
	}

	private Vector<NTT<PFE>> sampleBoundedVector(int dimension, byte[] seed, int nonce) {
		List<NTT<PFE>> result = new ArrayList<>();
		for (int i = 0; i < dimension; i++) {
			byte[] input = Arrays.copyOf(seed, seed.length + 2);
			input[seed.length] = (byte) (nonce + i);
			result.add(sampleBounded(hashH.stream(input)));
		}
		return new Vector<>(result);
	}

	private NTT<PFE> sampleUnbounded(InputStream stream) {
		List<PFE> coeffs = new ArrayList<>();
		Iterator<IntE> it = reinterpretStreamUnbounded(stream);
		while (coeffs.size() < degree) {
			coeffs.add(field.getInteger(it.next()));
		}
		return cylcotomic.fromList(coeffs);
	}

	public Matrix<NTT<PFE>> expandMatrix(byte[] seed) {
		List<List<NTT<PFE>>> result = new ArrayList<>();
		for (int i = 0; i < rows; i++) {
			List<NTT<PFE>> row = new ArrayList<>();
			for (int j = 0; j < cols; j++) {
				byte[] input = Arrays.copyOf(seed, seed.length + 2);
				input[seed.length] = (byte) j;
				input[seed.length + 1] = (byte) i;
				row.add(sampleUnbounded(hashG.stream(input)));
			}
			result.add(row);
		}
		return new Matrix<>(result);
	}

	private Vector<NTT<PFE>> expandMask(byte[] seed, int kappa) {
		List<NTT<PFE>> result = new ArrayList<>();
		for (int i = 0; i < cols; i++) {
			byte[] input = Arrays.copyOf(seed, seed.length + 2);
			input[seed.length] = (byte) ((kappa + i) % 256);
			input[seed.length + 1] = (byte) ((kappa + i) / 256);
			byte[] byteMask = hashH.evaluate(input, (log2Gamma1 + 1) * degree / 8);
			List<PFE> coeffs = new ArrayList<>();
			for (int j = 0; j < degree; j++) {
				int index = (log2Gamma1 + 1) * j;
				int endBitIndex = (log2Gamma1 + 1) * (j + 1);
				int value = 0;
				int power = 1;
				while (index < endBitIndex) {
					int mask = (1 << Math.min(8 - (index % 8), endBitIndex - index)) - 1;
					int bits = mask & (byteMask[index / 8] >> (index % 8));
					value += power * bits;
					power <<= 8 - (index % 8);
					index += 8 - (index % 8);
				}
				coeffs.add(field.getInteger(gamma1 - value));
			}
			result.add(cylcotomic.fromPolynomial(field.getUnivariatePolynomialRing().getPolynomial(coeffs)));
		}
		return new Vector<>(result);
	}

	private NTT<PFE> sampleInBall(byte[] seed) {
		InputStream stream = hashH.stream(seed);
		byte[] signs = new byte[8];
		try {
			stream.read(signs);
		} catch (IOException e) {
			throw new RuntimeException("Stream ended unexpectedly!", e);
		}
		List<PFE> result = new ArrayList<>();
		for (int i = 0; i < degree; i++) {
			result.add(field.zero());
		}
		for (int i = degree - tau; i < degree; i++) {
			int bitIndex = i - degree + tau;
			boolean sign = (signs[bitIndex / 8] & (1 << (bitIndex % 8))) != 0;
			int j;
			do {
				try {
					j = stream.read();
				} catch (IOException e) {
					throw new RuntimeException("Stream ended unexpectedly!", e);
				}
				if (j < 0) {
					throw new RuntimeException("Stream ended unexpectedly!");
				}
			} while (j > i);
			PFE tmp = result.get(j);
			result.set(j, sign ? field.getInteger(-1) : field.one());
			result.set(i, tmp);
		}
		return cylcotomic.fromPolynomial(field.getUnivariatePolynomialRing().getPolynomial(result));
	}

	private byte[] crhInput(byte[] prefix, Vector<IntE> t1) {
		IntPacker packer = new IntPacker(logQ - d);
		for (IntE c : t1.asList()) {
			packer.write(c.intValueExact());
		}
		return ByteArray.concatenate(prefix, packer.close()).array();
	}

	public byte[] serializePublicKey(DilithiumPublicKey publicKey) {
		return crhInput(publicKey.rho, publicKey.t1);
	}

	@Override
	public DilithiumPrivateKey createPrivateKey() {
		byte[] zeta = new byte[32];
		new SecureRandom().nextBytes(zeta);
		return createPrivateKeyFromSeed(zeta);
	}

	public DilithiumPrivateKey createPrivateKeyFromSeed(byte[] zeta) {
		byte[] input = Arrays.copyOf(zeta, zeta.length + 2);
		input[zeta.length] = (byte) rows;
		input[zeta.length + 1] = (byte) cols;
		byte[] expanded = hashH.evaluate(input, 128);
		byte[] rho = Arrays.copyOf(expanded, 32);
		byte[] sigma = Arrays.copyOfRange(expanded, 32, 3 * 32);
		byte[] k = Arrays.copyOfRange(expanded, 3 * 32, 4 * 32);
		Vector<NTT<PFE>> s = sampleBoundedVector(cols, sigma, 0);
		Vector<NTT<PFE>> e = sampleBoundedVector(rows, sigma, cols);
		return createPrivateKeyFromParameters(rho, k, s, e);
	}

	public DilithiumPrivateKey createPrivateKeyFromParameters(byte[] rho, byte[] k, Vector<NTT<PFE>> s,
			Vector<NTT<PFE>> e) {
		Matrix<NTT<PFE>> matrix = expandMatrix(rho);
		Vector<NTT<PFE>> t = freeModule.add(matrixModule.multiply(matrix, s), e);
		Pair<Vector<IntE>, Vector<NTT<PFE>>> power2Round = power2Round(t, d);
		byte[] tr = hashH.evaluate(crhInput(rho, power2Round.getFirst()), 64);
		return new DilithiumPrivateKey(new DilithiumPublicKey(rho, power2Round.getFirst()), k, tr, s, e,
				power2Round.getSecond());

	}

	@Override
	public DilithiumPublicKey createPublicKey(DilithiumPrivateKey privateKey) {
		return privateKey.getPublicKey();
	}

	@Override
	public DilithiumSignature sign(byte[] message, DilithiumPrivateKey privateKey) {
		Matrix<NTT<PFE>> matrix = expandMatrix(privateKey.getRho());
		byte[] input = ByteArray.concatenate(ByteArray.concatenate(privateKey.getTr(), new byte[2]).array(), message)
				.array();
		byte[] mu = hashH.evaluate(input, 64);
		int kappa = 0;
		Vector<NTT<PFE>> z = null;
		List<Boolean> hints = null;
		byte[] cSeed = null;
		ByteArray rhoPrimeInput = ByteArray.concatenate(ByteArray.concatenate(privateKey.k, new byte[32]).array(), mu);
		byte[] rhoPrime = hashH.evaluate(rhoPrimeInput.array(), 64);
		while (z == null) {
			Vector<NTT<PFE>> y = expandMask(rhoPrime, kappa);
			Vector<NTT<PFE>> w = matrixModule.multiply(matrix, y);
			Vector<IntE> w1 = highBits(w, 2 * gamma2);
			byte[] w1Packed = packVector(w1, w1bits);
			cSeed = hashH.evaluate(ByteArray.concatenate(mu, w1Packed).array(), lambda / 4);
			NTT<PFE> c = sampleInBall(cSeed);
			z = matrixModule.domain().add(y, matrixModule.domain().scalarMultiply(c, privateKey.getS()));
			Vector<NTT<PFE>> wcs2 = freeModule.subtract(w, freeModule.scalarMultiply(c, privateKey.getE()));
			Vector<NTT<PFE>> r0 = lowBits(wcs2, 2 * gamma2);
			if (infinityNorm(z) >= gamma1 - tau * eta || infinityNorm(r0) >= gamma2 - tau * eta) {
				z = null;
				cSeed = null;
			} else {
				Vector<NTT<PFE>> ct0 = freeModule.scalarMultiply(c, privateKey.getT0());
				hints = makeHints(freeModule.negative(ct0), freeModule.add(wcs2, ct0), 2 * gamma2);
				if (infinityNorm(ct0) >= gamma2 || oneNorm(hints) > omega) {
					z = null;
					hints = null;
					cSeed = null;
				}
			}
			kappa += cols;
		}
		return new DilithiumSignature(z, hints, cSeed);
	}

	// A method to create invalid signatures, that have a z value that is bigger
	// than the limit at the indicated index, but otherwise verify just fine
	// We do this by messing with the generation of y, to heavily put our thumb on the scale to make such a signature more likely, and then rejection sample
	public DilithiumSignature signViolateZ(byte[] message, DilithiumPrivateKey privateKey, int zIndexToViolate,
			boolean positive) {
		Matrix<NTT<PFE>> matrix = expandMatrix(privateKey.getRho());
		byte[] input = ByteArray.concatenate(ByteArray.concatenate(privateKey.getTr(), new byte[2]).array(), message)
				.array();
		byte[] mu = hashH.evaluate(input, 64);
		int kappa = 0;
		Vector<NTT<PFE>> z = null;
		List<Boolean> hints = null;
		byte[] cSeed = null;
		ByteArray rhoPrimeInput = ByteArray.concatenate(ByteArray.concatenate(privateKey.k, new byte[32]).array(), mu);
		byte[] rhoPrime = hashH.evaluate(rhoPrimeInput.array(), 64);
		int indexToViolate = zIndexToViolate / degree;
		PFE violationCoeff = field.getElement((positive ? 1 : -1) * (gamma1 - tau * eta / 2));
		while (z == null) {
			Vector<NTT<PFE>> y = expandMask(rhoPrime, kappa);
			UnivariatePolynomial<PFE> yInvNTT = cylcotomic.asPolynomial(y.get(indexToViolate + 1));
			List<PFE> coeffs = new ArrayList<>();
			for (int i = 0; i < degree; i++) {
				if (i == zIndexToViolate % degree) {
					coeffs.add(violationCoeff);
				} else {
					coeffs.add(yInvNTT.univariateCoefficient(i));
				}
			}
			NTT<PFE> violationYCoeff = cylcotomic.fromPolynomial(field.getUnivariatePolynomialRing().getPolynomial(coeffs));
			List<NTT<PFE>> modifiedYCoeffs = new ArrayList<>();
			for (int i = 0; i < cols; i++) {
				if (i == indexToViolate) {
					modifiedYCoeffs.add(violationYCoeff);
				} else {
					modifiedYCoeffs.add(y.get(i+1));
				}
			}
			y = new Vector<>(modifiedYCoeffs);
			Vector<NTT<PFE>> w = matrixModule.multiply(matrix, y);
			Vector<IntE> w1 = highBits(w, 2 * gamma2);
			byte[] w1Packed = packVector(w1, w1bits);
			cSeed = hashH.evaluate(ByteArray.concatenate(mu, w1Packed).array(), lambda / 4);
			NTT<PFE> c = sampleInBall(cSeed);
			z = matrixModule.domain().add(y, matrixModule.domain().scalarMultiply(c, privateKey.getS()));
			Vector<NTT<PFE>> wcs2 = freeModule.subtract(w, freeModule.scalarMultiply(c, privateKey.getE()));
			Vector<NTT<PFE>> r0 = lowBits(wcs2, 2 * gamma2);
			boolean reject = false;
			for (int i = 0; i < cols; i++) {
				UnivariatePolynomial<PFE> poly = cylcotomic.asPolynomial(z.get(i+1));
				if (reject) {
					break;
				}
				for (int j = 0; j < degree; j++) {
					if (i * degree + j == zIndexToViolate) {
						int value = Integers.z().centeredLift( poly.univariateCoefficient(j), q).intValueExact();
						if (value >= gamma1 || value < -gamma1) {
							reject = true;
							break;
						}
						if (positive && value < gamma1 - tau*eta) {
							reject = true;
							break;
						}
						if (!positive && value > -gamma1 + tau*eta) {
							reject = true;
							break;
						}
					} else {
						if (Integers.z().centeredLift( poly.univariateCoefficient(j), q).getValue().abs().intValueExact() >= gamma1 - tau*eta) {
							reject = true;
							break;
						}
					}
				}
			}
			if (reject || infinityNorm(r0) >= gamma2 - tau * eta) {
				z = null;
				cSeed = null;
			} else {
				Vector<NTT<PFE>> ct0 = freeModule.scalarMultiply(c, privateKey.getT0());
				hints = makeHints(freeModule.negative(ct0), freeModule.add(wcs2, ct0), 2 * gamma2);
				if (infinityNorm(ct0) >= gamma2 || oneNorm(hints) > omega) {
					z = null;
					hints = null;
					cSeed = null;
				}
			}
			kappa += cols;
		}
		return new DilithiumSignature(z, hints, cSeed);
	}

	public byte[] serializeSignature(DilithiumSignature signature) {
		Integers z = Integers.z();
		int size = 2 * gamma1 - 1;
		ByteArray result = new ByteArray(signature.cSeed);
		for (int i = 0; i < signature.z.dimension(); i++) {
			IntPacker packer = new IntPacker(MiscAlgorithms.bitLength(size));
			UnivariatePolynomial<PFE> coeff = cylcotomic.asPolynomial(signature.z.get(i + 1));
			for (int j = 0; j < degree; j++) {
				packer.write(gamma1 - z.centeredLift(coeff.univariateCoefficient(j), q).intValueExact());
			}
			result = ByteArray.concatenate(result.array(), packer.close());
		}
		byte[] hint = new byte[omega + rows];
		int index = 0;
		for (int i = 0; i < rows; i++) {
			for (int j = 0; j < degree; j++) {
				if (signature.hints.get(i * degree + j)) {
					hint[index] = (byte) j;
					index++;
				}
			}
			hint[i + omega] = (byte) index;
		}
		result = ByteArray.concatenate(result.array(), hint);
		return result.array();
	}

	@Override
	public boolean verify(byte[] message, DilithiumSignature signature, DilithiumPublicKey publicKey) {
		Matrix<NTT<PFE>> matrix = expandMatrix(publicKey.getRho());
		byte[] tr = hashH.evaluate(crhInput(publicKey.getRho(), publicKey.getT1()), 64);
		byte[] input = ByteArray.concatenate(ByteArray.concatenate(tr, new byte[2]).array(), message).array();
		byte[] mu = hashH.evaluate(input, 64);
		NTT<PFE> c = sampleInBall(signature.getCSeed());
		Vector<IntE> w1 = useHints(signature.getHints(),
				freeModule.subtract(matrixModule.multiply(matrix, signature.getZ()),
						freeModule.scalarMultiply(c, power2Inflate(publicKey.getT1(), d))),
				2 * gamma2);
		byte[] w1Packed = packVector(w1, w1bits);
		byte[] recoveredCSeed = hashH.evaluate(ByteArray.concatenate(mu, w1Packed).array(), lambda / 4);
		return infinityNorm(signature.getZ()) < gamma1 - tau * eta && Arrays.equals(recoveredCSeed, signature.cSeed)
				&& oneNorm(signature.getHints()) <= omega;
	}

	public boolean verifyNoCheck(byte[] message, DilithiumSignature signature, DilithiumPublicKey publicKey) {
		Matrix<NTT<PFE>> matrix = expandMatrix(publicKey.getRho());
		byte[] tr = hashH.evaluate(crhInput(publicKey.getRho(), publicKey.getT1()), 64);
		byte[] input = ByteArray.concatenate(ByteArray.concatenate(tr, new byte[2]).array(), message).array();
		byte[] mu = hashH.evaluate(input, 64);
		NTT<PFE> c = sampleInBall(signature.getCSeed());
		Vector<IntE> w1 = useHints(signature.getHints(),
				freeModule.subtract(matrixModule.multiply(matrix, signature.getZ()),
						freeModule.scalarMultiply(c, power2Inflate(publicKey.getT1(), d))),
				2 * gamma2);
		byte[] w1Packed = packVector(w1, w1bits);
		byte[] recoveredCSeed = hashH.evaluate(ByteArray.concatenate(mu, w1Packed).array(), lambda / 4);
		return Arrays.equals(recoveredCSeed, signature.cSeed) && oneNorm(signature.getHints()) <= omega;
	}
}
