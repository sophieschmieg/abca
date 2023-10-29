package cryptography;

import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.security.SecureRandom;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;

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
	private NTTRing<PFE, PrimeField> cylcotomic;
	private FreeModule<NTT<PFE>> freeModule;
	private MatrixModule<NTT<PFE>> matrixModule;
	private ExtendedOutputFunction xof;
	private ExtendedOutputFunction hash;
	private ExtendedOutputFunction crh;

	private final static int PRIME = (1 << 23) - (1 << 13) + 1;

	public final static Dilithium DILITHIUM2 = new Dilithium(4, 4, 2, 13, 39, 17, (PRIME - 1) / 88, 80, 4);
	public final static Dilithium DILITHIUM3 = new Dilithium(6, 5, 4, 13, 49, 19, (PRIME - 1) / 32, 55, 6);
	public final static Dilithium DILITHIUM5 = new Dilithium(8, 7, 2, 13, 60, 19, (PRIME - 1) / 32, 75, 6);

	public Dilithium(int k, int l, int eta, int d, int tau, int log2Gamma1, int gamma2, int omega, int w1bits) {
		this(PRIME, 256, 1753, k, l, eta, d, tau, log2Gamma1, gamma2, omega, w1bits, Sha3.SHAKE_128, Sha3.SHAKE_256,
				Sha3.SHAKE_256);
	}

	public Dilithium(int q, int degree, int r, int k, int l, int eta, int d, int tau, int log2Gamma1, int gamma2,
			int omega, int w1bits, ExtendedOutputFunction xof, ExtendedOutputFunction hash,
			ExtendedOutputFunction crh) {
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
		this.w1bits = w1bits;
		this.xof = xof;
		this.hash = hash;
		this.crh = crh;
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
			int logToWrite = log;
			while (logToWrite > 0) {
				int bits = Math.min(logToWrite, 8 - currentOffset);
				int mod = t % (1 << bits);
				current |= mod << currentOffset;
				currentOffset += bits;
				logToWrite -= bits;
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

	private Iterator<IntE> reinterpretStream(int log, InputStream stream) {
		return new Iterator<>() {
			private Integers z = Integers.z();
			private int bits = 0;
			private int leftOver = 0;
			private IntE next = null;

			private boolean getBit() throws IOException {
				if (bits == 0) {
					leftOver = stream.read();
					if (leftOver < 0) {
						throw new IOException("End of stream reached");
					}
					bits = 8;
				}
				boolean result = (leftOver & 0x01) != 0;
				bits--;
				leftOver >>= 1;
				return result;
			}

			@Override
			public boolean hasNext() {
				if (next != null) {
					return true;
				}
				next = z.zero();
				IntE power = z.one();
				try {
					for (int i = 0; i < log; i++) {
						if (getBit()) {
							next = z.add(next, power);
						}
						power = z.multiply(2, power);
					}
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

	private NTT<PFE> sampleUniform(IntE min, IntE max, InputStream stream) {
		Integers z = Integers.z();
		IntE range = z.add(z.subtract(max, min), z.one());
		int log = 0;
		IntE powerOfTwo = z.one();
		while (powerOfTwo.compareTo(range) < 0) {
			log++;
			powerOfTwo = z.multiply(2, powerOfTwo);
		}
		List<PFE> coeffs = new ArrayList<>();
		Iterator<IntE> it = reinterpretStream(log, stream);
		while (coeffs.size() < degree) {
			IntE candidate = it.next();
			if (candidate.compareTo(range) < 0) {
				coeffs.add(field.getInteger(z.add(candidate, min)));
			}
		}
		return cylcotomic.fromPolynomial(field.getUnivariatePolynomialRing().getPolynomial(coeffs));
	}

	private NTT<PFE> sampleUniform(IntE max, InputStream stream) {
		Integers z = Integers.z();
		return sampleUniform(z.negative(max), max, stream);
	}

	private Vector<NTT<PFE>> sampleUniformVector(IntE max, int dimension, byte[] seed, int nonce) {
		List<NTT<PFE>> result = new ArrayList<>();
		for (int i = 0; i < dimension; i++) {
			byte[] input = Arrays.copyOf(seed, seed.length + 1);
			input[seed.length] = (byte) (nonce + i);
			result.add(sampleUniform(max, hash.stream(input)));
		}
		return new Vector<>(result);
	}

	private NTT<PFE> sampleUniform(InputStream stream) {
		List<PFE> coeffs = new ArrayList<>();
		Iterator<IntE> it = reinterpretStream(logQ, stream);
		while (coeffs.size() < degree) {
			IntE candidate = it.next();
			if (candidate.compareTo(q) < 0) {
				coeffs.add(field.getInteger(candidate));
			}
		}
		return cylcotomic.fromList(coeffs);
	}

	private Matrix<NTT<PFE>> expandMatrix(byte[] seed) {
		List<List<NTT<PFE>>> result = new ArrayList<>();
		for (int i = 0; i < rows; i++) {
			List<NTT<PFE>> row = new ArrayList<>();
			for (int j = 0; j < cols; j++) {
				byte[] input = Arrays.copyOf(seed, seed.length + 2);
				input[seed.length] = (byte) j;
				input[seed.length + 1] = (byte) i;
				row.add(sampleUniform(xof.stream(input)));
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
			Iterator<IntE> mask = reinterpretStream(log2Gamma1, crh.stream(input));
			List<PFE> coeffs = new ArrayList<>();
			for (int j = 0; j < degree; j++) {
				coeffs.add(field.getInteger(mask.next()));
			}
			result.add(cylcotomic.fromPolynomial(field.getUnivariatePolynomialRing().getPolynomial(coeffs)));
		}
		return new Vector<>(result);
	}

	private NTT<PFE> sampleInBall(byte[] seed) {
		InputStream stream = crh.stream(seed);
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

	@Override
	public DilithiumPrivateKey createPrivateKey() {
		Integers z = Integers.z();
		byte[] zeta = new byte[32];
		new SecureRandom().nextBytes(zeta);
		byte[] expanded = hash.evaluate(zeta, 3 * 32);
		byte[] rho = Arrays.copyOf(expanded, 32);
		byte[] sigma = Arrays.copyOfRange(expanded, 32, 2 * 32);
		byte[] k = Arrays.copyOfRange(expanded, 2 * 32, 3 * 32);
		Matrix<NTT<PFE>> matrix = expandMatrix(rho);
		Vector<NTT<PFE>> s = sampleUniformVector(z.getInteger(eta), cols, sigma, 0);
		Vector<NTT<PFE>> e = sampleUniformVector(z.getInteger(eta), rows, sigma, cols);
		Vector<NTT<PFE>> t = freeModule.add(matrixModule.multiply(matrix, s), e);
		Pair<Vector<IntE>, Vector<NTT<PFE>>> power2Round = power2Round(t, d);
		byte[] tr = crh.evaluate(crhInput(rho, power2Round.getFirst()), 48);
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
		byte[] input = ByteArray.concatenate(privateKey.getTr(), message).array();
		byte[] mu = crh.evaluate(input, 48);
		int kappa = 0;
		Vector<NTT<PFE>> z = null;
		List<Boolean> hints = null;
		byte[] cSeed = null;
		byte[] rhoPrime = new byte[48];
		new SecureRandom().nextBytes(rhoPrime);
		while (z == null) {
			Vector<NTT<PFE>> y = expandMask(rhoPrime, kappa);
			Vector<NTT<PFE>> w = matrixModule.multiply(matrix, y);
			Vector<IntE> w1 = highBits(w, 2 * gamma2);
			byte[] w1Packed = packVector(w1, w1bits);
			cSeed = hash.evaluate(ByteArray.concatenate(mu, w1Packed).array(), 32);
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

	@Override
	public boolean verify(byte[] message, DilithiumSignature signature, DilithiumPublicKey publicKey) {
		Matrix<NTT<PFE>> matrix = expandMatrix(publicKey.getRho());
		byte[] tr = crh.evaluate(crhInput(publicKey.getRho(), publicKey.getT1()), 48);
		byte[] input = ByteArray.concatenate(tr, message).array();
		byte[] mu = crh.evaluate(input, 48);
		NTT<PFE> c = sampleInBall(signature.getCSeed());
		Vector<IntE> w1 = useHints(signature.getHints(),
				freeModule.subtract(matrixModule.multiply(matrix, signature.getZ()),
						freeModule.scalarMultiply(c, power2Inflate(publicKey.getT1(), d))),
				2 * gamma2);
		byte[] w1Packed = packVector(w1, w1bits);
		byte[] recoveredCSeed = hash.evaluate(ByteArray.concatenate(mu, w1Packed).array(), 32);
		return infinityNorm(signature.getZ()) < gamma1 - tau * eta && Arrays.equals(recoveredCSeed, signature.cSeed)
				&& oneNorm(signature.getHints()) <= omega;
	}
}
