package cryptography;

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
		private Vector<IntE> t0;

		public DilithiumPrivateKey(DilithiumPublicKey publicKey, byte[] k, byte[] tr, Vector<NTT<PFE>> s,
				Vector<NTT<PFE>> e, Vector<IntE> t0) {
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

		public Vector<IntE> getT0() {
			return t0;
		}

		@Override
		public int compareTo(DilithiumPrivateKey o) {
			return publicKey.compareTo(o.publicKey);
		}

	}

	public static class DilithiumSignature extends AbstractElement<DilithiumSignature> {
		private Vector<NTT<PFE>> z;
		private Vector<NTT<PFE>> h;
		private NTT<PFE> c;

		private DilithiumSignature(Vector<NTT<PFE>> z, Vector<NTT<PFE>> h, NTT<PFE> c) {
			this.z = z;
			this.h = h;
			this.c = c;
		}

		public Vector<NTT<PFE>> getZ() {
			return z;
		}

		public Vector<NTT<PFE>> getH() {
			return h;
		}

		public NTT<PFE> getC() {
			return c;
		}

		@Override
		public int compareTo(DilithiumSignature o) {
			int cmp = z.compareTo(o.z);
			if (cmp != 0) {
				return cmp;
			}
			cmp = h.compareTo(o.h);
			if (cmp != 0) {
				return cmp;
			}
			return c.compareTo(o.c);
		}

	}

	private IntE q;
	private int logQ;
	private PrimeField field;
	private int degree;
	private int rows;
	private int cols;
	private int eta;
	private NTTRing<PFE, PrimeField> cylcotomic;
	private MatrixModule<NTT<PFE>> matrixModule;
	private int d;
	private ExtendedOutputFunction xof;
	private ExtendedOutputFunction hash;
	private ExtendedOutputFunction crh;

	public Dilithium(int k, int l, int eta, int d) {
		this((1 << 23) - (1 << 13) + 1, 256, 1753, k, l, eta, d, Sha3.SHAKE_128, Sha3.SHAKE_256, Sha3.SHAKE_256);
	}

	public Dilithium(int q, int degree, int r, int k, int l, int eta, int d, ExtendedOutputFunction xof,
			ExtendedOutputFunction hash, ExtendedOutputFunction crh) {
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
		this.d = d;
		this.xof = xof;
		this.hash = hash;
		this.crh = crh;
	}

	private Vector<IntE> toIntVector(NTT<PFE> ntt) {
		List<PFE> asList = cylcotomic.asList(ntt);
		List<IntE> result = new ArrayList<>();
		for (PFE t : asList) {
			result.add(field.liftToInteger(t));
		}
		return new Vector<>(result);
	}

	private Vector<IntE> toIntVector(Vector<NTT<PFE>> ntt) {
		List<IntE> result = new ArrayList<>();
		for (int i = 0; i < ntt.dimension(); i++) {
			result.addAll(toIntVector(ntt.get(i + 1)).asList());
		}
		return new Vector<IntE>(result);
	}

	private NTT<PFE> toNTT(Vector<IntE> t) {
		return toNTT(t.asList());
	}

	private NTT<PFE> toNTT(List<IntE> t) {
		List<PFE> asCoeffList = new ArrayList<>();
		for (IntE c : t) {
			asCoeffList.add(field.getInteger(c));
		}
		return cylcotomic.fromPolynomial(field.getUnivariatePolynomialRing().getPolynomial(asCoeffList));
	}

	private Vector<NTT<PFE>> toNTTVector(Vector<IntE> t) {
		int dim = t.dimension() / degree;
		List<IntE> asList = t.asList();
		List<NTT<PFE>> result = new ArrayList<>();
		for (int i = 0; i < dim; i++) {
			result.add(toNTT(asList.subList(i * degree, (i + 1) * degree)));
		}
		return new Vector<>(result);
	}

	private Pair<Vector<IntE>, Vector<IntE>> power2Round(Vector<IntE> t, int exp) {
		List<IntE> result1 = new ArrayList<>();
		List<IntE> result2 = new ArrayList<>();
		for (int i = 0; i < t.dimension(); i++) {
			Pair<IntE, IntE> power2Round = power2Round(t.get(i + 1), exp);
			result1.add(power2Round.getFirst());
			result2.add(power2Round.getSecond());
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

	private Pair<IntE, IntE> decompose(IntE t, IntE alpha) {
		Integers z = Integers.z();
		t = z.remainder(t, q);
		IntE remainder = z.remainder(t, alpha);
		if (z.multiply(2, remainder).compareTo(alpha) > 0) {
			remainder = z.subtract(remainder, alpha);
		}
		if (z.subtract(t, remainder).equals(z.subtract(q, z.one()))) {
			return new Pair<>(z.zero(), z.subtract(remainder, z.one()));
		}
		IntE quotient = z.divideChecked(z.subtract(t, remainder), alpha);
		return new Pair<>(quotient, remainder);
	}

	private IntE highBits(IntE t, IntE alpha) {
		return decompose(t, alpha).getFirst();
	}

	private IntE lowBits(IntE t, IntE alpha) {
		return decompose(t, alpha).getSecond();
	}

	private boolean makeHint(IntE z, IntE r, IntE alpha) {
		IntE highR = highBits(r, alpha);
		IntE highV = highBits(Integers.z().add(z, r), alpha);
		return !highR.equals(highV);
	}

	private IntE useHint(boolean hint, IntE r, IntE alpha) {
		Integers z = Integers.z();
		IntE m = z.divideChecked(z.subtract(q, z.one()), alpha);
		Pair<IntE, IntE> decomposed = decompose(r, alpha);
		if (hint && decomposed.getSecond().compareTo(z.zero()) > 0) {
			return z.remainder(z.add(decomposed.getFirst(), z.one()), m);
		}
		if (hint && decomposed.getSecond().compareTo(z.zero()) <= 0) {
			return z.remainder(z.subtract(decomposed.getFirst(), z.one()), m);
		}
		return decomposed.getFirst();
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

	private byte[] crhInput(byte[] prefix, Vector<IntE> t1) {
		Integers z = Integers.z();
		byte[] input = new byte[MiscAlgorithms.DivRoundUp(8 * prefix.length + t1.dimension() * logQ, 8)];
		System.arraycopy(prefix, 0, input, 0, prefix.length);
		int byteIndex = prefix.length;
		int bitIndex = 0;
		for (int i = 0; i < t1.dimension(); i++) {
			IntE number = t1.get(i + 1);
			int bits = logQ;
			while (bits > 0) {
				int write = Math.min(8 - bitIndex, bits);
				int mod = 1 << write;
				int littleEndian = z.remainder(number, z.getInteger(mod)).intValueExact();
				input[byteIndex] = (byte) (littleEndian << bitIndex);
				bitIndex += write;
				if (bitIndex == 8) {
					bitIndex = 0;
					byteIndex++;
				}
				bits -= write;
				number = z.divideChecked(z.subtract(number, z.getInteger(littleEndian)), z.getInteger(mod));
			}
		}
		return input;
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
		Vector<NTT<PFE>> t = matrixModule.codomain().add(matrixModule.multiply(matrix, s), e);
		Vector<IntE> asIntVector = toIntVector(t);
		Pair<Vector<IntE>, Vector<IntE>> power2Round = power2Round(asIntVector, d);
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
		Vector<IntE> z = null;
		Vector<IntE> hints = null;
		byte[] rhoPrime = new byte[48];
		new SecureRandom().nextBytes(rhoPrime);
		while (z == null) {
			
		}
		return null;
	}

	@Override
	public boolean verify(byte[] message, DilithiumSignature signature, DilithiumPublicKey publicKey) {
		// TODO Auto-generated method stub
		return false;
	}
}
