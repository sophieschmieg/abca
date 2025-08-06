package cryptography;

import java.io.IOException;
import java.io.InputStream;
import java.security.SecureRandom;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import cryptography.Falcon.FalconPrivateKey;
import cryptography.Falcon.FalconPublicKey;
import cryptography.Falcon.FalconSignature;
import cryptography.interfaces.ExtendedOutputFunction;
import cryptography.interfaces.SignatureScheme;
import fields.finitefields.PrimeField;
import fields.finitefields.PrimeField.PFE;
import fields.floatingpoint.Complex;
import fields.floatingpoint.Complex.ComplexNumber;
import fields.floatingpoint.Reals;
import fields.floatingpoint.Reals.Real;
import fields.helper.AbstractElement;
import fields.helper.NTTRing;
import fields.helper.NTTRing.NTT;
import fields.integers.Integers;
import fields.integers.Integers.IntE;
import fields.integers.Rationals;
import fields.integers.Rationals.Fraction;
import fields.interfaces.Element;
import fields.interfaces.MathMap;
import fields.interfaces.Ring.ExtendedEuclideanResult;
import fields.interfaces.UnivariatePolynomial;
import fields.interfaces.UnivariatePolynomialRing;
import fields.numberfields.NumberField;
import fields.numberfields.NumberField.NFE;
import fields.vectors.FreeModule;
import fields.vectors.Matrix;
import fields.vectors.MatrixAlgebra;
import fields.vectors.Vector;
import util.FunctionMathMap;
import util.MiscAlgorithms;
import util.Pair;

public class Falcon implements SignatureScheme<FalconSignature, FalconPublicKey, FalconPrivateKey, Falcon> {
	
	public class FalconPublicKey extends AbstractElement<FalconPublicKey> {
		private NFE h;

		public FalconPublicKey(NFE h) {
			this.h = h;
		}

		@Override
		public int compareTo(FalconPublicKey o) {
			return h.compareTo(o.h);
		}

		public NFE getH() {
			return h;
		}

	}

	public class FalconPrivateKey extends AbstractElement<FalconPrivateKey> {
		private NFE f;
		private NFE g;
		private NFE capitalF;
		private NFE capitalG;
		private Matrix<NTT<ComplexNumber>> b;
		private FalconTree root;
		private FalconPublicKey publicKey;

		private FalconPrivateKey(NFE f, NFE g, NFE capitalF, NFE capitalG) {
			this.f = f;
			this.g = g;
			this.capitalF = capitalF;
			this.capitalG = capitalG;
			this.b = twoByTwo(asFFT(g), fftTower.get(log2degree).negative(asFFT(f)), asFFT(capitalG),
					fftTower.get(log2degree).negative(asFFT(capitalF)));
		}

		@Override
		public int compareTo(FalconPrivateKey o) {
			int cmp = f.compareTo(o.f);
			if (cmp != 0) {
				return cmp;
			}
			cmp = g.compareTo(o.g);
			if (cmp != 0) {
				return cmp;
			}
			cmp = capitalF.compareTo(o.capitalF);
			if (cmp != 0) {
				return cmp;
			}
			return capitalG.compareTo(o.capitalG);
		}

		public NFE getF() {
			return f;
		}

		public NFE getG() {
			return g;
		}

		public NFE getCapitalF() {
			return capitalF;
		}

		public NFE getCapitalG() {
			return capitalG;
		}

	}

	public class FalconSignature extends AbstractElement<FalconSignature> {
		private ByteArray random;
		private NFE signature;

		public FalconSignature(byte[] random, NFE signature) {
			this.random = new ByteArray(random);
			this.signature = signature;
		}

		@Override
		public int compareTo(FalconSignature o) {
			int cmp = random.compareTo(o.random);
			if (cmp != 0) {
				return cmp;
			}
			return signature.compareTo(o.signature);
		}

		public NFE getSignature() {
			return signature;
		}

	}

	private class FalconTree {
		private Real leafValue;
		private NTT<ComplexNumber> innerValue;
		private FalconTree left;
		private FalconTree right;
	}

	private int prime;
	private int degree;
	private int log2degree;
	private Complex c;
	private Reals r;
	private MathMap<Fraction, ComplexNumber> embeddingMap;
	private MathMap<ComplexNumber, Fraction> roundMap;
	private NTTRing<PFE, PrimeField> ntt;
	private List<NTTRing<ComplexNumber, Complex>> fftTower;
	private List<NumberField> nfTower;
	private List<List<ComplexNumber>> rootPowers;
	private SecureRandom random;
	private Real sigmaMin;
	private Real sigma;
	private Real sigmaMax;
	private Real limit;
	
	private static final int[][] rcdt = new int[][] { { 0x02, 0x18, 0x39, 0xac, 0xd3, 0x2e, 0xf4, 0xf7, 0xa3 },
			{ 0x82, 0xdb, 0x7d, 0x3f, 0x1f, 0x18, 0x2b, 0xd3, 0x54 },
			{ 0xff, 0xc1, 0x29, 0x48, 0x93, 0xd0, 0xcd, 0x7d, 0x22 },
			{ 0xe4, 0x4a, 0x99, 0xc7, 0x77, 0x43, 0x75, 0xd1, 0x0a },
			{ 0x6f, 0x1f, 0x3f, 0xf3, 0xae, 0x6c, 0x84, 0x95, 0x02 },
			{ 0x5f, 0xbd, 0x74, 0xed, 0x54, 0xc7, 0x4a, 0x77, 0x00 },
			{ 0xe4, 0x6a, 0x77, 0x2b, 0x54, 0xdd, 0x24, 0x10, 0x00 },
			{ 0xda, 0x63, 0xad, 0x65, 0xdc, 0xff, 0xa1, 0x01, 0x00 },
			{ 0x28, 0x64, 0x7b, 0x8a, 0xd8, 0x80, 0x1f, 0x00, 0x00 },
			{ 0x69, 0x0c, 0x04, 0xb2, 0xfd, 0xc3, 0x01, 0x00, 0x00 },
			{ 0xfb, 0x31, 0xd0, 0x24, 0xcf, 0x12, 0x00, 0x00, 0x00 },
			{ 0x1f, 0x09, 0x8b, 0x9f, 0x94, 0x00, 0x00, 0x00, 0x00 },
			{ 0x98, 0xa9, 0x5d, 0x66, 0x03, 0x00, 0x00, 0x00, 0x00 },
			{ 0xbb, 0x6e, 0xbf, 0x0e, 0x00, 0x00, 0x00, 0x00, 0x00 },
			{ 0x7e, 0x5d, 0x2f, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00 },
			{ 0x98, 0x70, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00 },
			{ 0xc6, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00 },
			{ 0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00 },
			{ 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00 } };

	public Falcon(int degree) {
		this(12289, degree, 1.277833697, 165.736617183, 1.8205, 34034726);
	}

	public Falcon(int prime, int log2degree, double sigmaMin, double sigma, double sigmaMax, double limit) {
		Integers z = Integers.z();
		this.c = Complex.c(64);
		this.r = c.getReals();
		this.random = new SecureRandom();
		this.prime = prime;
		this.degree = 1;
		this.log2degree = log2degree;
		this.sigmaMin = r.getDouble(sigmaMin);
		this.sigma = r.getDouble(sigma);
		this.sigmaMax = r.getDouble(sigmaMax);
		this.limit = r.getDouble(limit);
		this.fftTower = new ArrayList<>();
		this.rootPowers = new ArrayList<>();
		this.nfTower = new ArrayList<>();
		UnivariatePolynomialRing<IntE> zPolynomials = z.getUnivariatePolynomialRing();
		for (int i = 0; i <= log2degree; i++) {
			UnivariatePolynomial<IntE> mipo = zPolynomials.add(zPolynomials.getVarPower(degree), zPolynomials.one());
			this.nfTower.add(NumberField.getNumberFieldFromIntegerPolynomial(mipo));
			this.fftTower.add(new NTTRing<>(c, degree));
			this.rootPowers.add(new ArrayList<>());
			for (int j = 0; j < degree / 2; j++) {
				int power = 2 * NTTRing.bitReverse(j, i + 1) + 1;
				Vector<Real> cosineSine = r.cosineAndSine(r.divide(r.multiply(power, r.pi()), r.getInteger(degree)));
				this.rootPowers.get(i).add(c.fromVector(cosineSine));
			}
			if (i < log2degree) {
				degree *= 2;
			}

		}
		PFE root;
		PrimeField fp = PrimeField.getPrimeField(prime);
		do {
			root = fp.power(fp.getRandomElement(), (prime - 1) / (2 * degree));
		} while (!fp.power(root, degree).equals(fp.getInteger(-1)));
		this.ntt = new NTTRing<>(fp, degree, root);
		Rationals q = Rationals.q();
		this.embeddingMap = new FunctionMathMap<>((Fraction t) -> c.getFraction(t));
		this.roundMap = new FunctionMathMap<>((ComplexNumber t) -> q.getInteger(t.realPart().round()));
	}

	private <T extends Element<T>> Matrix<T> twoByTwo(T a11, T a12, T a21, T a22) {
		List<List<T>> asMatrix = new ArrayList<>();
		asMatrix.add(new ArrayList<>());
		asMatrix.add(new ArrayList<>());
		asMatrix.get(0).add(a11);
		asMatrix.get(0).add(a12);
		asMatrix.get(1).add(a21);
		asMatrix.get(1).add(a22);
		return new Matrix<>(asMatrix);
	}

	public NTT<ComplexNumber> asFFT(NFE t) {
		return asFFT(log2degree, t);
	}

	public NFE fromFFT(NTT<ComplexNumber> t) {
		return fromFFT(log2degree, t);
	}

	public NumberField getNumberField(int layer) {
		return nfTower.get(layer);
	}

	public NTTRing<ComplexNumber, Complex> getFFTRing(int layer) {
		return fftTower.get(layer);
	}

	public NTTRing<PFE, PrimeField> getNTTRing() {
		return ntt;
	}

	public int getPrime() {
		return prime;
	}
	
	public int getLog2Degree() {
		return log2degree;
	}
	
	public int getDegree() {
		return degree;
	}

	private NTT<ComplexNumber> asFFT(int layer, NFE t) {
		NTTRing<ComplexNumber, Complex> fft = fftTower.get(layer);
		return fft.fromPolynomial(c.getUnivariatePolynomialRing().getEmbedding(t.asPolynomial(), embeddingMap));
	}

	private NFE fromFFT(int layer, NTT<ComplexNumber> t) {
		Rationals q = Rationals.q();
		UnivariatePolynomial<Fraction> asPolynomial = q.getUnivariatePolynomialRing()
				.getEmbedding(fftTower.get(layer).asPolynomial(t), roundMap);
		return nfTower.get(layer).fromPolynomial(asPolynomial);
	}

	// t0 + t1*X, t0, t1 in Q[X^2]
	// t0 is even, t1 is odd coeffs of FFT
	private NTT<ComplexNumber> mergeFFT(int layer, NTT<ComplexNumber> t0, NTT<ComplexNumber> t1) {
		// t(zeta^(2*bitreverse(i)+1)) = t0(zeta^(2*(2*bitreverse(i)+1)) +
		// zeta^(2*bitreverse(i)+1)*t1(zeta^(2*(2*bitreverse(i)+1))
		// t0(zeta^2*bitreverse(i)+1)
		int step = 1 << (layer - 1);
		List<ComplexNumber> result = new ArrayList<>();
		for (int i = 0; i < step; i++) {
			ComplexNumber stepRoot = fftTower.get(layer).getRootPowers().get(i + step);
			ComplexNumber odd = c.multiply(stepRoot, t1.get(i));
			ComplexNumber even = t0.get(i);
			result.add(c.add(odd, even));
			result.add(c.subtract(even, odd));
		}
		return fftTower.get(layer).fromList(result);
	}

	private Vector<NTT<ComplexNumber>> splitFFT(int layer, NTT<ComplexNumber> t) {
		List<ComplexNumber> asList = t.values();
		List<ComplexNumber> t0 = new ArrayList<>();
		List<ComplexNumber> t1 = new ArrayList<>();
		int offset = 1 << (layer - 1);
		for (int i = 0; i < offset; i++) {
			ComplexNumber stepRoot = fftTower.get(layer).getInverseRootPowers().get(i + offset);
			ComplexNumber odd = asList.get(2 * i + 1);
			ComplexNumber even = asList.get(2 * i);
			t0.add(c.divide(c.add(odd, even), c.getInteger(2)));
			t1.add(c.divide(c.multiply(stepRoot, c.subtract(even, odd)), c.getInteger(2)));
		}
		return new Vector<>(fftTower.get(layer - 1).fromList(t0), fftTower.get(layer - 1).fromList(t1));
	}

	private int baseSampler() {
		byte[] u = new byte[9];
		random.nextBytes(u);
		int result = 0;
		// constant_time_compare? Never heard of her.
		outer: for (; result < 18; result++) {
			for (int i = 8; i >= 0; i--) {
				int value = u[i] & 255;
				if (value > rcdt[result][i]) {
					break outer;
				}
				if (value < rcdt[result][i]) {
					continue outer;
				}
			}
		}
		return result;
	}

	private IntE approxExp(Real x, Real ccs) {
		return r.multiply(r.power(r.getInteger(2), 63), ccs, r.exp(r.negative(x))).round();
	}

	private boolean berExp(Real x, Real ccs) {
		Integers z = Integers.z();
		IntE s = r.divide(x, r.log(r.getInteger(2))).roundDown();
		Real m = r.subtract(x, r.multiply(s, r.log(r.getInteger(2))));
		s = z.min(s, z.getInteger(63));
		IntE b = z.divide(z.subtract(z.multiply(2, approxExp(m, ccs)), z.one()), z.power(z.getInteger(2), s));
		int i = 64;
		int w;
		do {
			i -= 8;
			w = random.nextInt(256)
					- z.remainder(z.divide(b, z.power(z.getInteger(2), i)), z.getInteger(256)).intValueExact();
		} while (w == 0 && i >= 0);
		return w < 0;
	}

	private IntE samplerZ(Real mu, Real sigma) {
		Integers z = Integers.z();
		Real s = r.subtract(mu, r.getInteger(mu.roundDown()));
		Real ccs = r.divide(sigmaMin, sigma);
		while (true) {
			IntE z0 = z.getInteger(baseSampler());
			int b = random.nextBoolean() ? 1 : 0;
			IntE zx = z.add(z.getInteger(b), z.multiply(2 * b - 1, z0));
			Real x = r.subtract(r.divide(r.power(r.subtract(r.getInteger(zx), s), 2), r.multiply(2, sigma, sigma)),
					r.divide(r.power(r.getInteger(z0), 2), r.multiply(2, sigmaMax, sigmaMax)));
			if (berExp(x, ccs)) {
				return z.add(zx, mu.roundDown());
			}
		}
	}

	private Vector<NTT<ComplexNumber>> fastFourierSampler(int layer, NTT<ComplexNumber> t0, NTT<ComplexNumber> t1,
			FalconTree t) {
		if (layer == 0) {
			IntE z0 = samplerZ(t0.values().get(0).realPart(), t.leafValue);
			IntE z1 = samplerZ(t1.values().get(0).realPart(), t.leafValue);
			return new Vector<>(fftTower.get(layer).getInteger(z0), fftTower.get(layer).getInteger(z1));
		}
		NTTRing<ComplexNumber, Complex> complexRing = fftTower.get(layer);
		NTT<ComplexNumber> l = t.innerValue;
		FalconTree tree0 = t.left;
		FalconTree tree1 = t.right;
		Vector<NTT<ComplexNumber>> split1 = splitFFT(layer, t1);
		Vector<NTT<ComplexNumber>> splitZ1 = fastFourierSampler(layer - 1, split1.get(1), split1.get(2), tree1);
		NTT<ComplexNumber> z1 = mergeFFT(layer, splitZ1.get(1), splitZ1.get(2));
		NTT<ComplexNumber> t0Prime = complexRing.add(t0, complexRing.multiply(complexRing.subtract(t1, z1), l));
		Vector<NTT<ComplexNumber>> split0 = splitFFT(layer, t0Prime);
		Vector<NTT<ComplexNumber>> splitZ0 = fastFourierSampler(layer - 1, split0.get(1), split0.get(2), tree0);
		NTT<ComplexNumber> z0 = mergeFFT(layer, splitZ0.get(1), splitZ0.get(2));
		return new Vector<>(z0, z1);
	}

	private Pair<NFE, NFE> asLower(int layer, NFE t) {
		Rationals q = Rationals.q();
		UnivariatePolynomialRing<Fraction> poly = q.getUnivariatePolynomialRing();
		UnivariatePolynomial<Fraction> asPoly = t.asPolynomial();
		List<Fraction> even = new ArrayList<>();
		List<Fraction> odd = new ArrayList<>();
		for (int i = 0; i < (1 << layer); i++) {
			if (i % 2 == 0) {
				even.add(asPoly.univariateCoefficient(i));
			} else {
				odd.add(asPoly.univariateCoefficient(i));
			}
		}
		return new Pair<>(nfTower.get(layer - 1).fromPolynomial(poly.getPolynomial(even)),
				nfTower.get(layer - 1).fromPolynomial(poly.getPolynomial(odd)));
	}

	private NFE fromLower(int layer, Pair<NFE, NFE> t) {
		Rationals q = Rationals.q();
		UnivariatePolynomialRing<Fraction> poly = q.getUnivariatePolynomialRing();
		NFE even = nfTower.get(layer).fromPolynomial(
				poly.substitute(t.getFirst().asPolynomial(), Collections.singletonList(poly.getVarPower(2))));
		NFE odd = nfTower.get(layer).fromPolynomial(
				poly.substitute(t.getSecond().asPolynomial(), Collections.singletonList(poly.getVarPower(2))));
		return nfTower.get(layer).add(even, nfTower.get(layer).multiply(nfTower.get(layer).alpha(), odd));
	}

	private NFE norm(int layer, NFE t) {
		Pair<NFE, NFE> lower = asLower(layer, t);
		return nfTower.get(layer - 1).subtract(nfTower.get(layer - 1).power(lower.getFirst(), 2),
				nfTower.get(layer - 1).multiply(nfTower.get(layer - 1).alpha(), lower.getSecond(), lower.getSecond()));
	}

	private NFE conjugate(int layer, NFE t) {
		Pair<NFE, NFE> asLower = asLower(layer, t);
		return fromLower(layer, new Pair<>(asLower.getFirst(), nfTower.get(layer - 1).negative(asLower.getSecond())));
	}

	public Pair<NFE, NFE> reduce(int layer, NFE f, NFE g, NFE coeffF, NFE coeffG) {
		NFE k;
		Rationals q = Rationals.q();
		UnivariatePolynomialRing<Fraction> poly = q.getUnivariatePolynomialRing();
		Integers z = Integers.z();
		IntE min = z.zero();
		IntE max = z.zero();
		for (int i = 0; i < degree; i++) {
			min = z.min(f.asPolynomial().univariateCoefficient(i).asInteger(), min);
			min = z.min(g.asPolynomial().univariateCoefficient(i).asInteger(), min);
			max = z.max(f.asPolynomial().univariateCoefficient(i).asInteger(), max);
			max = z.max(g.asPolynomial().univariateCoefficient(i).asInteger(), max);
		}
		int mantissaBitLength = 53;
		NumberField nf = nfTower.get(layer);
		NTTRing<ComplexNumber, Complex> fft = fftTower.get(layer);
		int size = Math.max(mantissaBitLength, Math.max(min.getValue().bitLength(), max.getValue().bitLength()));
		NFE fAdj = nf.fromPolynomial(poly.getEmbedding(f.asPolynomial(), new FunctionMathMap<>(
				(Fraction t) -> q.getInteger(t.asInteger().getValue().shiftRight(size - mantissaBitLength)))));
		NFE gAdj = nf.fromPolynomial(poly.getEmbedding(g.asPolynomial(), new FunctionMathMap<>(
				(Fraction t) -> q.getInteger(t.asInteger().getValue().shiftRight(size - mantissaBitLength)))));
		NTT<ComplexNumber> fftF = asFFT(layer, fAdj);
		NTT<ComplexNumber> fftG = asFFT(layer, gAdj);
		NTT<ComplexNumber> denominator = fft.add(fft.multiply(fftF, adjoint(layer, fftF)),
				fft.multiply(fftG, adjoint(layer, fftG)));
		do {
			min = z.zero();
			max = z.zero();
			for (int i = 0; i < (1 << layer); i++) {
				min = z.min(coeffF.asPolynomial().univariateCoefficient(i).asInteger(), min);
				min = z.min(coeffG.asPolynomial().univariateCoefficient(i).asInteger(), min);
				max = z.max(coeffF.asPolynomial().univariateCoefficient(i).asInteger(), max);
				max = z.max(coeffG.asPolynomial().univariateCoefficient(i).asInteger(), max);
			}
			int coeffSize = Math.max(mantissaBitLength,
					Math.max(min.getValue().bitLength(), max.getValue().bitLength()));
			if (coeffSize < size) {
				break;
			}
			NFE coeffFAdj = nf.fromPolynomial(poly.getEmbedding(coeffF.asPolynomial(), new FunctionMathMap<>(
					(Fraction t) -> q.getInteger(t.asInteger().getValue().shiftRight(coeffSize - mantissaBitLength)))));
			NFE coeffGAdj = nf.fromPolynomial(poly.getEmbedding(coeffG.asPolynomial(), new FunctionMathMap<>(
					(Fraction t) -> q.getInteger(t.asInteger().getValue().shiftRight(coeffSize - mantissaBitLength)))));
			NTT<ComplexNumber> fftCoeffF = asFFT(layer, coeffFAdj);
			NTT<ComplexNumber> fftCoeffG = asFFT(layer, coeffGAdj);

			NTT<ComplexNumber> coeff = fft.divide(fft.add(fft.multiply(fftCoeffF, adjoint(layer, fftF)),
					fft.multiply(fftCoeffG, adjoint(layer, fftG))), denominator);
			k = fromFFT(layer, coeff);
			System.out
					.println("Current bit length: " + coeffSize + " target bit length: " + size + " multiplier: " + k);
			coeffF = nf.subtract(coeffF, nf.multiply(z.power(z.getInteger(2), coeffSize - size), k, f));
			coeffG = nf.subtract(coeffG, nf.multiply(z.power(z.getInteger(2), coeffSize - size), k, g));
		} while (!k.equals(nf.zero()));
		return new Pair<>(coeffF, coeffG);
	}

	public Pair<NFE, NFE> ntruSolve(int layer, NFE f, NFE g) {
		Integers z = Integers.z();
		NumberField nf = nfTower.get(layer);
		if (layer == 0) {
			IntE fc = f.asPolynomial().univariateCoefficient(0).asInteger();
			IntE gc = g.asPolynomial().univariateCoefficient(0).asInteger();
			ExtendedEuclideanResult<IntE> egcd = z.extendedEuclidean(fc, gc);
			IntE coeffF;
			IntE coeffG;
			if (egcd.getGcd().equals(z.getInteger(-1))) {
				coeffF = z.negative(egcd.getCoeff1());
				coeffG = egcd.getCoeff2();
			} else if (!egcd.getGcd().equals(z.one())) {
				return null;
			} else {
				coeffF = egcd.getCoeff1();
				coeffG = z.negative(egcd.getCoeff2());
			}
			return reduce(layer, f, g, nf.multiply(prime, nf.getInteger(coeffG)),
					nf.multiply(prime, nf.getInteger(coeffF)));
		}
		NFE normF = norm(layer, f);
		NFE normG = norm(layer, g);
		Pair<NFE, NFE> lower = ntruSolve(layer - 1, normF, normG);
		if (lower == null) {
			return null;
		}
		NFE coeffF = nf.multiply(fromLower(layer, new Pair<>(lower.getFirst(), nfTower.get(layer - 1).zero())),
				conjugate(layer, g));
		NFE coeffG = nf.multiply(fromLower(layer, new Pair<>(lower.getSecond(), nfTower.get(layer - 1).zero())),
				conjugate(layer, f));
		return reduce(layer, f, g, coeffF, coeffG);
	}

	public NFE adjoint(int layer, NFE t) {
		if (layer == 0) {
			return t;
		}
		Rationals q = Rationals.q();
		UnivariatePolynomialRing<Fraction> rationalPolynomialRing = q.getUnivariatePolynomialRing();
		List<Fraction> coeffs = new ArrayList<>();
		coeffs.add(t.asPolynomial().univariateCoefficient(0));
		for (int i = 1; i < (1 << layer); i++) {
			coeffs.add(q.negative(t.asPolynomial().univariateCoefficient(degree - i)));
		}
		return nfTower.get(layer).fromPolynomial(rationalPolynomialRing.getPolynomial(coeffs));
	}

	public NTT<ComplexNumber> adjoint(int layer, NTT<ComplexNumber> t) {
		List<ComplexNumber> result = new ArrayList<>();
		for (ComplexNumber coeff : t.values()) {
			result.add(c.conjugate(coeff));
		}
		return fftTower.get(layer).fromList(result);
	}

	private FalconPrivateKey ntruKeyGen() {
		Integers z = Integers.z();
		UnivariatePolynomialRing<IntE> intPolynomialRing = z.getUnivariatePolynomialRing();
		Rationals q = Rationals.q();
		UnivariatePolynomialRing<Fraction> fractionsPolynomialRing = q.getUnivariatePolynomialRing();
		Real sigma = r.multiply(1.17, r.positiveSqrt(r.divide(r.getInteger(prime), r.getInteger(8192))));
		NumberField nf = nfTower.get(log2degree);
		while (true) {
			List<IntE> fCoeffs = new ArrayList<>();
			List<IntE> gCoeffs = new ArrayList<>();
			for (int i = 0; i < degree; i++) {
				IntE fCoeff = z.zero();
				IntE gCoeff = z.zero();
				for (int j = 0; j < 4096 / degree; j++) {
					fCoeff = z.add(fCoeff, samplerZ(r.zero(), sigma));
					gCoeff = z.add(gCoeff, samplerZ(r.zero(), sigma));
				}
				fCoeffs.add(fCoeff);
				gCoeffs.add(gCoeff);
			}
			NFE f = nf.fromPolynomial(fractionsPolynomialRing.getEmbedding(intPolynomialRing.getPolynomial(fCoeffs),
					q.getEmbeddingMap()));
			NFE g = nf.fromPolynomial(fractionsPolynomialRing.getEmbedding(intPolynomialRing.getPolynomial(gCoeffs),
					q.getEmbeddingMap()));
			NTT<PFE> reducedF = reduceAndTransform(f);
			NTT<PFE> reducedG = reduceAndTransform(g);
			if (!ntt.isUnit(reducedF) || !ntt.isUnit(reducedG)) {
				continue;
			}
			NTTRing<ComplexNumber, Complex> ring = fftTower.get(log2degree);
			NTT<ComplexNumber> denominator = ring.add(ring.multiply(asFFT(f), asFFT(adjoint(log2degree, f))),
					ring.multiply(asFFT(g), asFFT(adjoint(log2degree, g))));
			Real gamma = r.max(normSqr(log2degree, new Vector<>(f, nf.negative(g))),
					normSqr(new Vector<>(ring.divide(asFFT(nf.multiply(prime, adjoint(log2degree, f))), denominator),
							ring.divide(asFFT(nf.multiply(prime, adjoint(log2degree, g))), denominator))));
			if (gamma.compareTo(r.multiply(1.17 * 1.17, r.getInteger(prime))) > 0) {
				continue;
			}
			Pair<NFE, NFE> capitalFG = ntruSolve(log2degree, f, g);
			if (capitalFG == null) {
				continue;
			}
			return new FalconPrivateKey(f, g, capitalFG.getFirst(), capitalFG.getSecond());
		}
	}

	private Real innerProduct(int layer, NFE t1, NFE t2) {
		Rationals q = Rationals.q();
		Fraction result = q.zero();
		for (int i = 0; i < 1 << layer; i++) {
			result = q.add(
					q.multiply(t1.asPolynomial().univariateCoefficient(i), t2.asPolynomial().univariateCoefficient(i)),
					result);
		}
		return r.getEmbedding(result);
	}

	private Real innerProduct(int layer, Vector<NFE> t1, Vector<NFE> t2) {
		Real result = r.zero();
		for (int i = 0; i < t1.dimension(); i++) {
			result = r.add(innerProduct(layer, t1.get(i + 1), t2.get(i + 1)), result);
		}
		return result;
	}

	private ComplexNumber innerProduct(NTT<ComplexNumber> t1, NTT<ComplexNumber> t2) {
		ComplexNumber result = c.zero();
		for (int i = 0; i < degree; i++) {
			result = c.add(c.multiply(t1.get(i), c.conjugate(t2.get(i))), result);
		}
		return result;
	}

	private ComplexNumber innerProduct(Vector<NTT<ComplexNumber>> t1, Vector<NTT<ComplexNumber>> t2) {
		ComplexNumber result = c.zero();
		for (int i = 0; i < t1.dimension(); i++) {
			result = c.add(innerProduct(t1.get(i + 1), t2.get(i + 1)), result);
		}
		return c.divide(result, c.getInteger(degree));
	}

	public Real normSqr(int layer, Vector<NFE> t) {
		return innerProduct(layer, t, t);
	}

	public Real normSqr(Vector<NTT<ComplexNumber>> t) {
		return c.asReal(innerProduct(t, t));
	}

	private Pair<Matrix<NTT<ComplexNumber>>, Matrix<NTT<ComplexNumber>>> ldl(int layer, Matrix<NTT<ComplexNumber>> g) {
		NTTRing<ComplexNumber, Complex> ring = fftTower.get(layer);
		NTT<ComplexNumber> d00 = g.entry(1, 1);
		NTT<ComplexNumber> l10 = ring.divideChecked(g.entry(2, 1), d00);
		NTT<ComplexNumber> d11 = ring.subtract(g.entry(2, 2), ring.multiply(l10, adjoint(layer, l10), d00));
		return new Pair<>(twoByTwo(ring.one(), ring.zero(), l10, ring.one()),
				twoByTwo(d00, ring.zero(), ring.zero(), d11));
	}

	private FalconTree ffldl(int layer, Matrix<NTT<ComplexNumber>> g) {
		Pair<Matrix<NTT<ComplexNumber>>, Matrix<NTT<ComplexNumber>>> ldl = ldl(layer, g);
		FalconTree node = new FalconTree();
		node.innerValue = ldl.getFirst().entry(2, 1);
		if (layer == 1) {
			node.left = new FalconTree();
			node.left.leafValue = r.multiply(sigma, r.inverseSqrt(ldl.getSecond().entry(1, 1).get(0).realPart()));
			node.right = new FalconTree();
			node.right.leafValue = r.multiply(sigma, r.inverseSqrt(ldl.getSecond().entry(2, 2).get(0).realPart()));
			return node;
		}
		Vector<NTT<ComplexNumber>> d0 = splitFFT(layer, ldl.getSecond().entry(1, 1));
		Vector<NTT<ComplexNumber>> d1 = splitFFT(layer, ldl.getSecond().entry(2, 2));
		Matrix<NTT<ComplexNumber>> g0 = twoByTwo(d0.get(1), d0.get(2), adjoint(layer - 1, d0.get(2)), d0.get(1));
		Matrix<NTT<ComplexNumber>> g1 = twoByTwo(d1.get(1), d1.get(2), adjoint(layer - 1, d1.get(2)), d1.get(1));
		node.left = ffldl(layer - 1, g0);
		node.right = ffldl(layer - 1, g1);
		return node;
	}

	public NTT<PFE> reduceAndTransform(NFE t) {
		PrimeField fp = PrimeField.getPrimeField(prime);
		UnivariatePolynomialRing<PFE> ring = fp.getUnivariatePolynomialRing();
		return ntt.fromPolynomial(
				ring.getEmbedding(t.asPolynomial(), new FunctionMathMap<>((Fraction s) -> fp.getFraction(s))));
	}

	public NFE lift(NTT<PFE> t) {
		Rationals q = Rationals.q();
		Integers z = Integers.z();
		return nfTower.get(log2degree).fromPolynomial(q.getUnivariatePolynomialRing().getEmbedding(ntt.asPolynomial(t),
				new FunctionMathMap<PFE, Fraction>((PFE s) -> q.getInteger(z.centeredLift(s, prime)))));
	}

	private NFE hashToPoint(ByteArray seed) {
		Rationals q = Rationals.q();
		List<Fraction> coeffs = new ArrayList<>();
		int k = MiscAlgorithms.DivRoundDown(1 << 16, prime);
		ExtendedOutputFunction xof = Sha3.SHAKE_256;
		InputStream in = xof.stream(seed.array());
		int i = 0;
		try {
			while (i < degree) {
				int coeff = 0;
				int read = in.read();
				if (read < 0) {
					throw new IOException();
				}
				coeff += read;
				read = in.read();
				if (read < 0) {
					throw new IOException();
				}
				coeff += 256 * read;
				if (coeff < k * prime) {
					coeffs.add(q.getInteger(coeff));
					i++;
				}
			}
		} catch (IOException e) {
			throw new RuntimeException("XOF ran out of bytes!", e);
		}
		return nfTower.get(log2degree).fromPolynomial(q.getUnivariatePolynomialRing().getPolynomial(coeffs));
	}

	@Override
	public FalconPublicKey createPublicKey(FalconPrivateKey privateKey) {
		return privateKey.publicKey;
	}

	public FalconPrivateKey createPrivateKeyFromGiven(NFE f, NFE g) {
		Pair<NFE, NFE> capitalFG = ntruSolve(log2degree, f, g);
		if (capitalFG == null) {
			throw new ArithmeticException("Could not solve!");
		}
		FalconPrivateKey key = new FalconPrivateKey(f, g, capitalFG.getFirst(), capitalFG.getSecond());
		Matrix<NTT<ComplexNumber>> b = key.b;
		FreeModule<NTT<ComplexNumber>> module = new FreeModule<>(fftTower.get(log2degree), 2);
		MatrixAlgebra<NTT<ComplexNumber>> algebra = module.matrixAlgebra();
		Matrix<NTT<ComplexNumber>> fftBt = algebra.transpose(
				Matrix.mapMatrix(new FunctionMathMap<>((NTT<ComplexNumber> s) -> adjoint(log2degree, s)), b));
		Matrix<NTT<ComplexNumber>> gMatrix = algebra.multiply(b, fftBt);
		key.root = ffldl(log2degree, gMatrix);
		key.publicKey = new FalconPublicKey(
				lift(ntt.divideChecked(reduceAndTransform(key.g), reduceAndTransform(key.f))));
		return key;
		
	}
	
	@Override
	public FalconPrivateKey createPrivateKey() {
		FalconPrivateKey key = ntruKeyGen();
		Matrix<NTT<ComplexNumber>> b = key.b;
		FreeModule<NTT<ComplexNumber>> module = new FreeModule<>(fftTower.get(log2degree), 2);
		MatrixAlgebra<NTT<ComplexNumber>> algebra = module.matrixAlgebra();
		Matrix<NTT<ComplexNumber>> fftBt = algebra.transpose(
				Matrix.mapMatrix(new FunctionMathMap<>((NTT<ComplexNumber> s) -> adjoint(log2degree, s)), b));
		Matrix<NTT<ComplexNumber>> g = algebra.multiply(b, fftBt);
		key.root = ffldl(log2degree, g);
		key.publicKey = new FalconPublicKey(
				lift(ntt.divideChecked(reduceAndTransform(key.g), reduceAndTransform(key.f))));
		return key;
	}

	public Pair<ByteArray, NTT<ComplexNumber>> hash(byte[] message) {
		byte[] random = new byte[40];
		this.random.nextBytes(random);
		NTT<ComplexNumber> c = asFFT(hashToPoint(ByteArray.concatenate(random, message)));
		return new Pair<>(new ByteArray(random), c);
	}

	public FalconSignature prehashSign(FalconPrivateKey privateKey, byte[] random, NTT<ComplexNumber> c) {
		NTTRing<ComplexNumber, Complex> ring = fftTower.get(log2degree);
		NTT<ComplexNumber> qInv = ring.inverse(ring.getInteger(prime));
		FreeModule<NTT<ComplexNumber>> module = new FreeModule<>(ring, 2);
		Vector<NTT<ComplexNumber>> t = new Vector<>(ring.multiply(-1, qInv, c, asFFT(privateKey.capitalF)),
				ring.multiply(qInv, c, asFFT(privateKey.f)));
		Vector<NTT<ComplexNumber>> s;
		do {
			Vector<NTT<ComplexNumber>> intSample = fastFourierSampler(log2degree, t.get(1), t.get(2), privateKey.root);
			s = module.matrixAlgebra().multiply(module.matrixAlgebra().transpose(privateKey.b),
					module.subtract(t, intSample));
		} while (normSqr(s).compareTo(limit) > 0);
		NFE signature = fromFFT(s.get(2));
		return new FalconSignature(random, signature);
	}

	@Override
	public FalconSignature sign(byte[] message, FalconPrivateKey privateKey) {
		Pair<ByteArray, NTT<ComplexNumber>> hashed = hash(message);
		return prehashSign(privateKey, hashed.getFirst().array(), hashed.getSecond());
	}

	@Override
	public boolean verify(byte[] message, FalconSignature signature, FalconPublicKey publicKey) {
		NFE c = hashToPoint(ByteArray.concatenate(signature.random.array(), message));
		NFE s1 = lift(ntt.subtract(reduceAndTransform(c),
				ntt.multiply(reduceAndTransform(publicKey.h), reduceAndTransform(signature.signature))));
		return innerProduct(log2degree, s1, signature.signature).compareTo(limit) <= 0;
	}
}
