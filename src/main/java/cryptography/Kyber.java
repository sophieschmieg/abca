package cryptography;

import java.io.IOException;
import java.io.InputStream;
import java.security.SecureRandom;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import cryptography.interfaces.ExtendedOutputFunction;
import cryptography.interfaces.HashFunction;
import cryptography.interfaces.KemScheme;
import cryptography.interfaces.PseudoRandomFunction;
import fields.finitefields.FiniteField;
import fields.finitefields.FiniteField.FFE;
import fields.finitefields.PrimeField;
import fields.finitefields.PrimeField.PFE;
import fields.helper.NTTRing;
import fields.helper.NTTRing.NTT;
import fields.integers.Integers;
import fields.integers.Integers.IntE;
import fields.integers.Rationals;
import fields.interfaces.Polynomial;
import fields.interfaces.UnivariatePolynomial;
import fields.interfaces.UnivariatePolynomialRing;
import fields.vectors.FreeModule;
import fields.vectors.Matrix;
import fields.vectors.MatrixAlgebra;
import fields.vectors.Vector;
import util.FunctionMathMap;
import util.MiscAlgorithms;
import util.Pair;

public class Kyber implements KemScheme<VariableLengthKey, ByteArray, ByteArray, ByteArray, Kyber> {
	private Integers z;
	private Rationals q;
	private int prime;
	private int log2Prime;
	private PrimeField field;
	private FiniteField extension;
	private int degree;
	private int degreeInBytes;
	private int k;
	private int eta1;
	private int eta2;
	private int du;
	private int dv;
	private HashFunction g;
	private HashFunction h;
	private PseudoRandomFunction prf;
	private ExtendedOutputFunction xof;
	private ExtendedOutputFunction kdf;
	private NTTRing<FFE, FiniteField> cyclotomic;
	private FreeModule<NTT<FFE>> freeModule;
	private MatrixAlgebra<NTT<FFE>> matrixAlgebra;
	private UnivariatePolynomialRing<PFE> polynomialRing;
	private UnivariatePolynomialRing<IntE> intPolynomialRing;
	private FreeModule<Polynomial<PFE>> freePolynomialModule;

	public Kyber(int k, int eta1, int eta2, int du, int dv, HashFunction g, HashFunction h, PseudoRandomFunction prf,
			ExtendedOutputFunction xof, ExtendedOutputFunction kdf) {
		this(256, 3329, 17, k, eta1, eta2, du, dv, g, h, prf, xof, kdf);
	}

	public Kyber(int degree, int prime, int zeta, int k, int eta1, int eta2, int du, int dv, HashFunction g,
			HashFunction h, PseudoRandomFunction prf, ExtendedOutputFunction xof, ExtendedOutputFunction kdf) {
		this.z = Integers.z();
		this.q = Rationals.q();
		this.prime = prime;
		this.log2Prime = 0;
		int powerOfTwo = 1;
		while (powerOfTwo < prime) {
			powerOfTwo *= 2;
			this.log2Prime++;
		}
		this.degree = degree;
		if (degree % 8 != 0) {
			throw new ArithmeticException("Degree needs to be a multiple of 8!");
		}
		this.degreeInBytes = degree / 8;
		this.k = k;
		this.eta1 = eta1;
		this.eta2 = eta2;
		this.du = du;
		this.dv = dv;
		this.g = g;
		this.h = h;
		this.prf = prf;
		this.xof = xof;
		this.kdf = kdf;
		this.field = PrimeField.getPrimeField(prime);
		this.intPolynomialRing = z.getUnivariatePolynomialRing();
		this.polynomialRing = field.getUnivariatePolynomialRing();
		this.extension = FiniteField.getFiniteField(
				polynomialRing.subtract(polynomialRing.getVarPower(2), polynomialRing.getInteger(zeta)), field);
		this.cyclotomic = new NTTRing<>(extension, degree, extension.alpha());
		this.freeModule = new FreeModule<>(cyclotomic, k);
		this.matrixAlgebra = freeModule.matrixAlgebra();
		this.freePolynomialModule = new FreeModule<>(polynomialRing, k);
	}

	private NTT<FFE> fromPolynomial(UnivariatePolynomial<PFE> polynomial) {
		return cyclotomic.fromPolynomial(
				extension.getUnivariatePolynomialRing().getEmbedding(polynomial, extension.getEmbeddingMap()));
	}

	private Vector<NTT<FFE>> fromPolynomialVector(Vector<Polynomial<PFE>> t) {
		return Vector.mapVector(
				new FunctionMathMap<>((Polynomial<PFE> s) -> fromPolynomial(polynomialRing.toUnivariate(s))), t);
	}

	private NTT<FFE> nttFromList(List<Integer> list) {
		List<FFE> values = new ArrayList<>();
		for (int i = 0; i < degree / 2; i++) {
			values.add(extension.fromVector(
					new Vector<>(field.getInteger(list.get(2 * i)), field.getInteger(list.get(2 * i + 1)))));
			values.add(extension.fromVector(
					new Vector<>(field.getInteger(list.get(2 * i)), field.getInteger(-list.get(2 * i + 1)))));
		}
		return cyclotomic.fromList(values);
	}

	private List<Integer> nttAsList(NTT<FFE> t) {
		List<FFE> asList = cyclotomic.asList(t);
		List<Integer> values = new ArrayList<>();
		for (int i = 0; i < degree / 2; i++) {
			Vector<PFE> asVector = extension.asVector(asList.get(2 * i));
			values.add(asVector.get(1).getValue().intValueExact());
			values.add(asVector.get(2).getValue().intValueExact());
			Vector<PFE> check = extension.asVector(asList.get(2 * i + 1));
			if (!check.get(1).equals(asVector.get(1)) || !field.negative(check.get(2)).equals(asVector.get(2))) {
				throw new ArithmeticException("List did not contain complements!");
			}
		}
		return values;
	}

	private UnivariatePolynomial<IntE> intPolynomialFromList(List<Integer> list) {
		List<IntE> coeffs = new ArrayList<>();
		for (int c : list) {
			coeffs.add(z.getInteger(c));
		}
		return z.getUnivariatePolynomialRing().getPolynomial(coeffs);
	}

	private List<Integer> intPolynomialAsList(UnivariatePolynomial<IntE> t) {
		List<Integer> coeffs = new ArrayList<>();
		for (int i = 0; i < degree; i++) {
			coeffs.add(t.univariateCoefficient(i).intValueExact());
		}
		return coeffs;
	}

	private UnivariatePolynomial<PFE> asPolynomial(NTT<FFE> t) {
		return polynomialRing.getEmbedding(cyclotomic.asPolynomial(t), extension.asBaseFieldElementMap());
	}

	private Vector<Polynomial<PFE>> asPolynomialVector(Vector<NTT<FFE>> t) {
		return Vector.mapVector(new FunctionMathMap<>((NTT<FFE> s) -> asPolynomial(s)), t);
	}

	private NTT<FFE> parse(InputStream in) throws IOException {
		List<Integer> result = new ArrayList<>();
		while (result.size() < degree) {
			if (log2Prime == 12) {
				int bPrev = in.read();
				int b = in.read();
				int bNext = in.read();
				if (bPrev < 0 || b < 0 || bNext < 0) {
					throw new IOException("Stream ended unexpectedly!");
				}
				int d1 = bPrev + 256 * (b % 16);
				int d2 = b / 16 + 16 * bNext;
				if (d1 < prime) {
					result.add(d1);
				}
				if (d2 < prime && result.size() < degree) {
					result.add(d2);
				}
			} else {
				int numBytes = MiscAlgorithms.DivRoundUp(log2Prime, 8);
				int integer = 0;
				for (int i = 0; i < numBytes; i++) {
					int b = in.read();
					if (b < 0) {
						throw new IOException("Stream ended unexpectedly!");
					}
					integer *= 256;
					integer += b;
				}
				integer &= (1 << log2Prime) - 1;
				if (integer < prime) {
					result.add(integer);
				}
			}
		}
		return nttFromList(result);
	}

	private UnivariatePolynomial<PFE> centeredBinominalDistribution(int eta, InputStream in) {
		byte[] dataArray = new byte[eta * 2 * degreeInBytes];
		for (int i = 0; i < dataArray.length; i++) {
			try {
				int read = in.read();
				if (read < 0) {
					throw new ArithmeticException("Stream ended!");
				}
				dataArray[i] = (byte) read;
			} catch (IOException e) {
				throw new RuntimeException("input stream had read error!", e);
			}
		}
		ByteArray data = new ByteArray(dataArray);
		List<PFE> result = new ArrayList<>();
		for (int i = 0; i < degree; i++) {
			PFE f = field.zero();
			for (int j = 0; j < eta; j++) {
				if (data.bitAt(2 * i * eta + j)) {
					f = field.add(f, field.one());
				}
				if (data.bitAt((2 * i + 1) * eta + j)) {
					f = field.subtract(f, field.one());
				}
			}
			result.add(f);
		}
		return polynomialRing.getPolynomial(result);
	}

	private ByteArray encode(List<Integer> t, int l) {
		byte[] encoded = new byte[degreeInBytes * l];
		int index = 0;
		int bitIndex = 0;
		for (int i = 0; i < degree; i++) {
			int value = t.get(i);
			if (value < 0 || value >= 1 << l) {
				throw new ArithmeticException("Cannot encode polynomial value + " + value + " in " + l + " bits!");
			}
			for (int j = 0; j < l; j++) {
				if ((value & (1 << j)) != 0) {
					encoded[index] |= 1 << bitIndex;
				}
				if (bitIndex == 7) {
					bitIndex = 0;
					index++;
				} else {
					bitIndex++;
				}
			}
		}
		return new ByteArray(encoded);
	}

	private List<Integer> decode(byte[] data) {
		if (data.length % degreeInBytes != 0) {
			throw new ArithmeticException("Not an encoded polynomial!");
		}
		int l = data.length / degreeInBytes;
		ByteArray t = new ByteArray(data);
		List<Integer> coeffs = new ArrayList<>();
		for (int i = 0; i < degree; i++) {
			int f = 0;
			for (int j = 0; j < l; j++) {
				if (t.bitAt(i * l + j)) {
					f |= 1 << j;
				}
			}
			coeffs.add(f);
		}
		return coeffs;
	}

	private ByteArray encodeNTT(NTT<FFE> t, int l) {
		return encode(nttAsList(t), l);
	}

	private ByteArray encodeIntPolynomial(UnivariatePolynomial<IntE> t, int l) {
		return encode(intPolynomialAsList(t), l);
	}

	private ByteArray encodeNTTVector(Vector<NTT<FFE>> t, int l) {
		byte[] result = new byte[degreeInBytes * l * t.dimension()];
		for (int i = 0; i < t.dimension(); i++) {
			System.arraycopy(encodeNTT(t.get(i + 1), l).array(), 0, result, degreeInBytes * l * i, degreeInBytes * l);
		}
		return new ByteArray(result);
	}

	private ByteArray encodeIntPolynomialVector(Vector<Polynomial<IntE>> t, int l) {
		byte[] result = new byte[degreeInBytes * l * t.dimension()];
		for (int i = 0; i < t.dimension(); i++) {
			System.arraycopy(encodeIntPolynomial(intPolynomialRing.toUnivariate(t.get(i + 1)), l).array(), 0, result,
					degreeInBytes * l * i, degreeInBytes * l);
		}
		return new ByteArray(result);
	}

	private NTT<FFE> decodeNTT(byte[] t) {
		return nttFromList(decode(t));
	}

	private UnivariatePolynomial<IntE> decodeIntPolynomial(byte[] t) {
		return intPolynomialFromList(decode(t));
	}

	private Vector<NTT<FFE>> decodeNTTVector(byte[] t, int l) {
		List<NTT<FFE>> elements = new ArrayList<>();
		int size = t.length / (degreeInBytes * l);
		for (int i = 0; i < size; i++) {
			byte[] copy = Arrays.copyOfRange(t, degreeInBytes * l * i, degreeInBytes * l * (i + 1));
			elements.add(decodeNTT(copy));
		}
		return new Vector<>(elements);
	}

	private Vector<Polynomial<IntE>> decodeIntPolynomialVector(byte[] t, int l) {
		List<Polynomial<IntE>> elements = new ArrayList<>();
		int size = t.length / (degreeInBytes * l);
		for (int i = 0; i < size; i++) {
			byte[] copy = Arrays.copyOfRange(t, degreeInBytes * l * i, degreeInBytes * l * (i + 1));
			elements.add(decodeIntPolynomial(copy));
		}
		return new Vector<>(elements);
	}

	private IntE compress(PFE t, int exp) {
		IntE power = z.power(z.getInteger(2), exp);
		return z.remainder(q.getFraction(z.multiply(power, field.liftToInteger(t)), z.getInteger(prime)).round(),
				power);
	}

	private PFE decompress(IntE t, int exp) {
		IntE power = z.power(z.getInteger(2), exp);
		return field.reduce(q.getFraction(z.multiply(prime, t), power).round());
	}

	private UnivariatePolynomial<IntE> compress(Polynomial<PFE> t, int exp) {
		return intPolynomialRing.getEmbedding(t, new FunctionMathMap<>((PFE s) -> compress(s, exp)));
	}

	private UnivariatePolynomial<PFE> decompress(Polynomial<IntE> t, int exp) {
		return polynomialRing.getEmbedding(t, new FunctionMathMap<>((IntE s) -> decompress(s, exp)));
	}

	private Vector<Polynomial<IntE>> compress(Vector<Polynomial<PFE>> t, int exp) {
		return Vector.mapVector(new FunctionMathMap<>((Polynomial<PFE> s) -> compress(s, exp)), t);
	}

	private Vector<Polynomial<PFE>> decompress(Vector<Polynomial<IntE>> t, int exp) {
		return Vector.mapVector(new FunctionMathMap<>((Polynomial<IntE> s) -> decompress(s, exp)), t);
	}

	private ByteArray encrypt(byte[] publicKey, byte[] random, byte[] message) {
		byte[] rho = Arrays.copyOfRange(publicKey, log2Prime * degreeInBytes * k, publicKey.length);
		Vector<NTT<FFE>> t = decodeNTTVector(Arrays.copyOf(publicKey, log2Prime * degreeInBytes * k), log2Prime);
		byte N = 0;
		List<List<NTT<FFE>>> matrixElements = new ArrayList<>();
		for (int i = 0; i < k; i++) {
			List<NTT<FFE>> row = new ArrayList<>();
			for (int j = 0; j < k; j++) {
				try {
					row.add(parse(xof.stream(
							ByteArray.concatenate(new ByteArray(rho), new ByteArray(new byte[] { (byte) i, (byte) j }))
									.array())));
				} catch (IOException e) {
					throw new RuntimeException("Extended output function has IOError!", e);
				}
			}
			matrixElements.add(row);
		}
		Matrix<NTT<FFE>> matrix = new Matrix<>(matrixElements);
		List<NTT<FFE>> rVectorElements = new ArrayList<>();
		for (int i = 0; i < k; i++) {
			rVectorElements
					.add(fromPolynomial(centeredBinominalDistribution(eta1, prf.stream(random, new byte[] { N }))));
			N++;
		}
		Vector<NTT<FFE>> r = new Vector<>(rVectorElements);
		List<Polynomial<PFE>> e1VectorElements = new ArrayList<>();
		for (int i = 0; i < k; i++) {
			e1VectorElements.add(centeredBinominalDistribution(eta2, prf.stream(random, new byte[] { N })));
			N++;
		}
		Vector<Polynomial<PFE>> e1 = new Vector<>(e1VectorElements);
		Polynomial<PFE> e2 = centeredBinominalDistribution(eta2, prf.stream(random, new byte[] { N }));
		Vector<Polynomial<PFE>> u = freePolynomialModule.add(asPolynomialVector(matrixAlgebra.multiply(matrix, r)), e1);
		UnivariatePolynomial<PFE> v = polynomialRing.toUnivariate(polynomialRing
				.add(asPolynomial(freeModule.innerProduct(t, r)), e2, decompress(decodeIntPolynomial(message), 1)));
		ByteArray c1 = encodeIntPolynomialVector(compress(u, du), du);
		ByteArray c2 = encodeIntPolynomial(compress(v, dv), dv);
		return ByteArray.concatenate(c1, c2);
	}

	private ByteArray decrypt(byte[] privateKey, byte[] ciphertext) {
		Vector<Polynomial<PFE>> u = decompress(
				decodeIntPolynomialVector(Arrays.copyOf(ciphertext, degreeInBytes * du * k), du), du);
		UnivariatePolynomial<PFE> v = decompress(
				decodeIntPolynomial(Arrays.copyOfRange(ciphertext, degreeInBytes * du * k, ciphertext.length)), dv);
		Vector<NTT<FFE>> s = decodeNTTVector(Arrays.copyOf(privateKey, log2Prime * degreeInBytes * k), log2Prime);
		return encodeIntPolynomial(compress(
				polynomialRing.subtract(v, asPolynomial(freeModule.innerProduct(s, fromPolynomialVector(u)))), 1), 1);
	}

	@Override
	public ByteArray createPublicKey(ByteArray privateKey) {
		return new ByteArray(Arrays.copyOfRange(privateKey.array(), log2Prime * degreeInBytes * k,
				log2Prime * degreeInBytes * k + log2Prime * degreeInBytes * k + degreeInBytes));
	}

	@Override
	public ByteArray createPrivateKey() {
		SecureRandom random = new SecureRandom();
		byte[] d = new byte[degreeInBytes];
		random.nextBytes(d);
		byte[] rhoAndSigma = g.evaluate(d);
		byte[] rho = Arrays.copyOf(rhoAndSigma, degreeInBytes);
		byte[] sigma = Arrays.copyOfRange(rhoAndSigma, degreeInBytes, 2 * degreeInBytes);
		byte N = 0;
		List<List<NTT<FFE>>> matrixElements = new ArrayList<>();
		for (int i = 0; i < k; i++) {
			List<NTT<FFE>> row = new ArrayList<>();
			for (int j = 0; j < k; j++) {
				try {
					row.add(parse(xof.stream(
							ByteArray.concatenate(new ByteArray(rho), new ByteArray(new byte[] { (byte) j, (byte) i }))
									.array())));
				} catch (IOException e) {
					throw new RuntimeException("Extended output function has IOError!", e);
				}
			}
			matrixElements.add(row);
		}
		Matrix<NTT<FFE>> matrix = new Matrix<>(matrixElements);
		List<NTT<FFE>> sVectorElements = new ArrayList<>();
		for (int i = 0; i < k; i++) {
			sVectorElements
					.add(fromPolynomial(centeredBinominalDistribution(eta1, prf.stream(sigma, new byte[] { N }))));
			N++;
		}
		Vector<NTT<FFE>> s = new Vector<>(sVectorElements);
		List<NTT<FFE>> eVectorElements = new ArrayList<>();
		for (int i = 0; i < k; i++) {
			eVectorElements
					.add(fromPolynomial(centeredBinominalDistribution(eta1, prf.stream(sigma, new byte[] { N }))));
			N++;
		}
		Vector<NTT<FFE>> e = new Vector<>(eVectorElements);
		Vector<NTT<FFE>> t = freeModule.add(matrixAlgebra.multiply(matrix, s), e);
		ByteArray publicKey = ByteArray.concatenate(encodeNTTVector(t, log2Prime), new ByteArray(rho));
		byte[] z = new byte[degreeInBytes];
		random.nextBytes(z);
		ByteArray privateKey = ByteArray.concatenate(ByteArray.concatenate(encodeNTTVector(s, log2Prime), publicKey),
				ByteArray.concatenate(new ByteArray(h.evaluate(publicKey.array())), new ByteArray(z)));
		return privateKey;
	}

	@Override
	public Pair<VariableLengthKey, ByteArray> encapsulate(ByteArray publicKey) {
		SecureRandom random = new SecureRandom();
		byte[] seed = new byte[degreeInBytes];
		random.nextBytes(seed);
		seed = Arrays.copyOf(h.evaluate(seed), degreeInBytes);
		byte[] preKeyAndRandomness = g.evaluate(ByteArray.concatenate(seed, h.evaluate(publicKey.array())).array());
		byte[] preKey = Arrays.copyOf(preKeyAndRandomness, degreeInBytes);
		byte[] randomness = Arrays.copyOfRange(preKeyAndRandomness, degreeInBytes, 2 * degreeInBytes);
		ByteArray ciphertext = encrypt(publicKey.array(), randomness, seed);
		return new Pair<>(
				new VariableLengthKey(ByteArray.concatenate(preKey, h.evaluate(ciphertext.array())).array(), kdf),
				ciphertext);
	}

	@Override
	public VariableLengthKey decapsulate(ByteArray kem, ByteArray privateKey) {
		byte[] publicKey = Arrays.copyOfRange(privateKey.array(), log2Prime * degreeInBytes * k,
				log2Prime * degreeInBytes * k + log2Prime * degreeInBytes * k + degreeInBytes);
		byte[] hash = Arrays.copyOfRange(privateKey.array(), 2 * log2Prime * degreeInBytes * k + degreeInBytes,
				2 * log2Prime * degreeInBytes * k + 2 * degreeInBytes);
		byte[] z = Arrays.copyOfRange(privateKey.array(), 2 * log2Prime * degreeInBytes * k + 2 * degreeInBytes,
				2 * log2Prime * degreeInBytes * k + 3 * degreeInBytes);
		byte[] recoveredSeed = decrypt(privateKey.array(), kem.array()).array();
		byte[] preKeyAndRandomness = g.evaluate(ByteArray.concatenate(recoveredSeed, hash).array());
		byte[] preKey = Arrays.copyOf(preKeyAndRandomness, degreeInBytes);
		byte[] randomness = Arrays.copyOfRange(preKeyAndRandomness, degreeInBytes, 2 * degreeInBytes);
		ByteArray ciphertext = encrypt(publicKey, randomness, recoveredSeed);
		if (ciphertext.equals(kem)) {
			return new VariableLengthKey(ByteArray.concatenate(preKey, h.evaluate(ciphertext.array())).array(), kdf);
		}
		return new VariableLengthKey(ByteArray.concatenate(z, h.evaluate(ciphertext.array())).array(), kdf);
	}
}
