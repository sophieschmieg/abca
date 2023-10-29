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
import fields.helper.FieldEmbedding;
import fields.helper.NTTRing;
import fields.helper.NTTRing.NTT;
import fields.helper.ProductRing;
import fields.helper.ProductRing.ProductElement;
import fields.integers.Integers;
import fields.integers.Integers.IntE;
import fields.integers.Rationals;
import fields.integers.Rationals.Fraction;
import fields.interfaces.Element;
import fields.interfaces.MathMap;
import fields.interfaces.Ring.ExtendedEuclideanResult;
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

	}

	public class FalconPrivateKey extends AbstractElement<FalconPrivateKey> {
		private NFE f;
		private NFE g;
		private NFE capitalF;
		private NFE capitalG;
		private Matrix<ProductElement<ComplexNumber>> b;
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

	}

	private class FalconTree {
		private Real leafValue;
		private ProductElement<ComplexNumber> innerValue;
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
	private NumberField nf;
	private NTTRing<PFE, PrimeField> ntt;
	private NTTRing<ComplexNumber, Complex> fft;
	private List<ProductRing<ComplexNumber, Complex>> fftTower;
	private List<ProductRing<IntE, Integers>> zFftTower;
	private List<FieldEmbedding<Fraction, NFE, NumberField>> nfTower;
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
		this.zFftTower = new ArrayList<>();
		this.nfTower = new ArrayList<>();
		this.rootPowers = new ArrayList<>();
		NumberField base = NumberField.getNumberField();
		for (int i = 0; i < log2degree; i++) {
			NFE alpha = base.degree() == 1 ? base.one() : base.alpha();
			UnivariatePolynomialRing<NFE> polynomialRing = base.getUnivariatePolynomialRing();
			this.nfTower.add(base.getEmbeddedExtension(
					polynomialRing.add(polynomialRing.getVarPower(2), polynomialRing.getEmbedding(alpha))));
			this.fftTower.add(ProductRing.power(c, degree));
			this.zFftTower.add(ProductRing.power(z, degree));
			this.rootPowers.add(new ArrayList<>());
			for (int j = 0; j < degree / 2; j++) {
				int power = 2 * NTTRing.bitReverse(j, i + 1) + 1;
				Vector<Real> cosineSine = r.cosineAndSine(r.divide(r.multiply(power, r.pi()), r.getInteger(degree)));
				this.rootPowers.get(i).add(c.fromVector(cosineSine));
			}
			degree *= 2;

		}
		this.fft = new NTTRing<>(c, degree);
		PFE root;
		PrimeField fp = PrimeField.getPrimeField(prime);
		do {
			root = fp.power(fp.getRandomElement(), (prime - 1) / (2 * degree));
		} while (!fp.power(root, degree).equals(fp.getInteger(-1)));
		this.ntt = new NTTRing<>(fp, degree, root);
		UnivariatePolynomialRing<IntE> intPolynomials = z.getUnivariatePolynomialRing();
		this.nf = NumberField.getNumberFieldFromIntegerPolynomial(
				intPolynomials.add(intPolynomials.getVarPower(degree), intPolynomials.one()));
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

	private ProductElement<ComplexNumber> asFFT(NFE t) {
		return fft.asProductElement(
				fft.fromPolynomial(c.getUnivariatePolynomialRing().getEmbedding(t.asPolynomial(), embeddingMap)));
	}

	private NFE fromFFT(NTT<ComplexNumber> t) {
		Rationals q = Rationals.q();
		return nf.fromPolynomial(q.getUnivariatePolynomialRing().getEmbedding(fft.asPolynomial(t), roundMap));
	}

//	private MathMap<NFE, ProductElement<ComplexNumber>> asFFT() {
//		return new FunctionMathMap<>((NFE t) -> asFFT(t));
//	}
//
//	private MathMap<NTT<ComplexNumber>, NFE> fromFFT() {
//		return new FunctionMathMap<>((NTT<ComplexNumber> t) -> fromFFT(t));
//	}

	private ProductElement<IntE> mergeFFTInt(int layer, ProductElement<IntE> t0, ProductElement<IntE> t1) {
		ProductRing<ComplexNumber, Complex> complexRing = fftTower.get(layer);
		ProductElement<ComplexNumber> z0 = complexRing.mapElement(t0,
				new FunctionMathMap<>((IntE s) -> c.getInteger(s)));
		ProductElement<ComplexNumber> z1 = complexRing.mapElement(t1,
				new FunctionMathMap<>((IntE s) -> c.getInteger(s)));
		ProductRing<IntE, Integers> intRing = zFftTower.get(layer + 1);
		return intRing.mapElement(mergeFFT(layer, z0, z1),
				new FunctionMathMap<>((ComplexNumber s) -> s.realPart().round()));
	}

	// t0 + t1*X, t0, t1 in Q[X^2]
	// t0 is even, t1 is odd coeffs of FFT
	private ProductElement<ComplexNumber> mergeFFT(int layer, ProductElement<ComplexNumber> t0,
			ProductElement<ComplexNumber> t1) {
		// t(zeta^(2*bitreverse(i)+1)) = t0(zeta^(2*(2*bitreverse(i)+1)) +
		// zeta^(2*bitreverse(i)+1)*t1(zeta^(2*(2*bitreverse(i)+1))
		// t0(zeta^2*bitreverse(i)+1)
		int step = 1 << layer;
		List<ComplexNumber> result = new ArrayList<>();
		for (int i = 0; i < step; i++) {
			ComplexNumber stepRoot = fft.getRootPowers().get(i + step);
			ComplexNumber odd = c.multiply(stepRoot, t1.get(i));
			ComplexNumber even = t0.get(i);
			result.add(c.add(odd, even));
			result.add(c.subtract(even, odd));
		}
		return fftTower.get(layer).getElement(result);
	}

	private Vector<ProductElement<ComplexNumber>> splitFFT(int layer, ProductElement<ComplexNumber> t) {
		List<ComplexNumber> asList = t.values();
		List<ComplexNumber> t0 = new ArrayList<>();
		List<ComplexNumber> t1 = new ArrayList<>();
		int offset = 1 << layer;
		for (int i = 0; i < offset; i++) {
			ComplexNumber stepRoot = fft.getInverseRootPowers().get(i + offset);
			ComplexNumber odd = asList.get(2 * i + 1);
			ComplexNumber even = asList.get(2 * i);
			t0.add(c.divide(c.add(odd, even), c.getInteger(2)));
			t1.add(c.divide(c.multiply(stepRoot, c.subtract(even, odd)), c.getInteger(2)));
		}
		return new Vector<>(fftTower.get(layer - 1).getElement(t0), fftTower.get(layer - 1).getElement(t1));
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
			IntE zx = z.add(z.getInteger(b), z.multiply(2 * (b - 1), z0));
			Real x = r.subtract(r.divide(r.power(r.subtract(r.getInteger(zx), s), 2), r.multiply(2, sigma, sigma)),
					r.divide(r.power(r.getInteger(z0), 2), r.multiply(2, sigmaMax, sigmaMax)));
			if (berExp(x, ccs)) {
				return z.add(zx, mu.roundDown());
			}
		}
	}

	private Pair<ProductElement<IntE>, ProductElement<IntE>> fastFourierSampler(int layer,
			ProductElement<ComplexNumber> t0, ProductElement<ComplexNumber> t1, FalconTree t) {
		if (layer == 0) {
			IntE z0 = samplerZ(t0.values().get(0).realPart(), t.leafValue);
			IntE z1 = samplerZ(t1.values().get(0).realPart(), t.leafValue);
			return new Pair<>(zFftTower.get(layer).getElement(Collections.singletonList(z0)),
					zFftTower.get(layer).getElement(Collections.singletonList(z1)));
		}
		ProductRing<ComplexNumber, Complex> complexRing = fftTower.get(layer);
		ProductElement<ComplexNumber> l = t.innerValue;
		FalconTree tree0 = t.left;
		FalconTree tree1 = t.right;
		Vector<ProductElement<ComplexNumber>> split1 = splitFFT(layer - 1, t1);
		Pair<ProductElement<IntE>, ProductElement<IntE>> splitZ1 = fastFourierSampler(layer - 1, split1.get(1),
				split1.get(2), tree1);
		ProductElement<IntE> z1 = mergeFFTInt(layer - 1, splitZ1.getFirst(), splitZ1.getSecond());
		ProductElement<ComplexNumber> t0Prime = complexRing.add(t0, complexRing.multiply(complexRing.subtract(t1,
				complexRing.mapElement(z1, new FunctionMathMap<>((IntE s) -> c.getInteger(s)))), l));
		Vector<ProductElement<ComplexNumber>> split0 = splitFFT(layer - 1, t0Prime);
		Pair<ProductElement<IntE>, ProductElement<IntE>> splitZ0 = fastFourierSampler(layer - 1, split0.get(1),
				split0.get(2), tree0);
		ProductElement<IntE> z0 = mergeFFTInt(layer - 1, splitZ0.getFirst(), splitZ0.getSecond());
		return new Pair<>(z0, z1);
	}

	private Vector<ProductElement<ComplexNumber>> asComplexVector(Pair<ProductElement<IntE>, ProductElement<IntE>> t) {
		List<ComplexNumber> first = new ArrayList<>();
		for (IntE s : t.getFirst().values()) {
			first.add(c.getInteger(s));
		}
		List<ComplexNumber> second = new ArrayList<>();
		for (IntE s : t.getSecond().values()) {
			second.add(c.getInteger(s));
		}
		return new Vector<>(fftTower.get(log2degree).getElement(first), fftTower.get(log2degree).getElement(second));
	}

	private Pair<NFE, NFE> ntruSolve(int layer, NFE f, NFE g) {
		Integers z = Integers.z();
		if (layer == 0) {
			FieldEmbedding<Fraction, NFE, NumberField> nf = nfTower.get(layer);
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
			return new Pair<>(nf.getEmbeddedField().multiply(prime, nf.getEmbeddedField().getInteger(coeffG)),
					nf.getEmbeddedField().multiply(prime, nf.getEmbeddedField().getInteger(coeffF)));
		}
		FieldEmbedding<Fraction, NFE, NumberField> nf = nfTower.get(layer - 1);
		NFE normF = nf.norm(f);
		NFE normG = nf.norm(g);
		Pair<NFE, NFE> lower = ntruSolve(layer - 1, normF, normG);
		if (lower == null) {
			return null;
		}
		UnivariatePolynomialRing<Fraction> qPolynomialRing = Rationals.q().getUnivariatePolynomialRing();
		NFE coeffF = nf.getField().fromPolynomial(qPolynomialRing.multiply(
				qPolynomialRing.substitute(lower.getFirst().asPolynomial(),
						Collections.singletonList(qPolynomialRing.getVarPower(2))),
				qPolynomialRing.substitute(g.asPolynomial(),
						Collections.singletonList(qPolynomialRing.getEmbedding(Rationals.q().getInteger(-1), 1)))));
		NFE coeffG = nf.getField().fromPolynomial(qPolynomialRing.multiply(
				qPolynomialRing.substitute(lower.getSecond().asPolynomial(),
						Collections.singletonList(qPolynomialRing.getVarPower(2))),
				qPolynomialRing.substitute(f.asPolynomial(),
						Collections.singletonList(qPolynomialRing.getEmbedding(Rationals.q().getInteger(-1), 1)))));
		NFE k;
		do {
			NFE coeff = nf.getField()
					.divide(nf.getField().add(nf.getField().multiply(coeffF, adjoint(layer, f)),
							nf.getField().multiply(coeffG, adjoint(layer, g))),
							nf.getField().add(nf.getField().multiply(f, adjoint(layer, f)),
									nf.getField().multiply(g, adjoint(layer, g))));
			k = nf.getField().fromPolynomial(qPolynomialRing.getEmbedding(coeff.asPolynomial(),
					new FunctionMathMap<>((Fraction t) -> Rationals.q().getInteger(t.round()))));
			coeffF = nf.getField().subtract(coeffF, nf.getField().multiply(k, f));
			coeffG = nf.getField().subtract(coeffG, nf.getField().multiply(k, g));
		} while (!k.equals(nf.getField().zero()));
		return new Pair<>(coeffF, coeffG);
	}

	private NFE adjoint(int layer, NFE t) {
		if (layer == 0) {
			return t;
		}
		NumberField nf = nfTower.get(layer - 1).getField();
		Rationals q = Rationals.q();
		UnivariatePolynomialRing<Fraction> rationalPolynomialRing = q.getUnivariatePolynomialRing();
		List<Fraction> coeffs = new ArrayList<>();
		coeffs.add(t.asPolynomial().univariateCoefficient(0));
		int degree = nf.degree();
		for (int i = 1; i < degree; i++) {
			coeffs.add(q.negative(t.asPolynomial().univariateCoefficient(degree - i)));
		}
		return nf.fromPolynomial(rationalPolynomialRing.getPolynomial(coeffs));
	}

	private ProductElement<ComplexNumber> adjoint(int layer, ProductElement<ComplexNumber> t) {
		List<ComplexNumber> result = new ArrayList<>();
		for (ComplexNumber coeff : t.values()) {
			result.add(c.conjugate(coeff));
		}
		return fftTower.get(layer).getElement(result);
	}

	private FalconPrivateKey ntruKeyGen() {
		Integers z = Integers.z();
		UnivariatePolynomialRing<IntE> intPolynomialRing = z.getUnivariatePolynomialRing();
		Rationals q = Rationals.q();
		UnivariatePolynomialRing<Fraction> fractionsPolynomialRing = q.getUnivariatePolynomialRing();
		Real sigma = r.multiply(1.17, r.positiveSqrt(r.divide(r.getInteger(prime), r.getInteger(8192))));
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
			ProductElement<ComplexNumber> denominator = asFFT(
					nf.add(nf.multiply(f, adjoint(log2degree, f)), nf.multiply(g, adjoint(log2degree, g))));
			ProductRing<ComplexNumber, Complex> ring = fftTower.get(log2degree);
			Real gamma = r.max(norm(log2degree, new Vector<>(f, nf.negative(g))),
					norm(new Vector<>(ring.divide(ring.multiply(prime, adjoint(log2degree, asFFT(f))), denominator),
							ring.divide(ring.multiply(prime, adjoint(log2degree, asFFT(g))), denominator))));
			if (gamma.compareTo(r.multiply(1.17, r.positiveSqrt(r.getInteger(prime)))) > 0) {
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

	private ComplexNumber innerProduct(ProductElement<ComplexNumber> t1, ProductElement<ComplexNumber> t2) {
		ComplexNumber result = c.zero();
		for (int i = 0; i < degree; i++) {
			result = c.add(c.multiply(t1.get(i), c.conjugate(t2.get(i))), result);
		}
		return result;
	}

	private ComplexNumber innerProduct(Vector<ProductElement<ComplexNumber>> t1,
			Vector<ProductElement<ComplexNumber>> t2) {
		ComplexNumber result = c.zero();
		for (int i = 0; i < t1.dimension(); i++) {
			result = c.add(innerProduct(t1.get(i + 1), t2.get(i + 1)), result);
		}
		return result;
	}

	private Real norm(int layer, Vector<NFE> t) {
		return r.positiveSqrt(innerProduct(layer, t, t));
	}

	private Real norm(Vector<ProductElement<ComplexNumber>> t) {
		return c.asReal(innerProduct(t, t));
	}

	private Pair<Matrix<ProductElement<ComplexNumber>>, Matrix<ProductElement<ComplexNumber>>> ldl(int layer,
			Matrix<ProductElement<ComplexNumber>> g) {
		ProductRing<ComplexNumber, Complex> ring = fftTower.get(layer);
		ProductElement<ComplexNumber> d00 = g.entry(1, 1);
		ProductElement<ComplexNumber> l10 = ring.divideChecked(g.entry(2, 1), d00);
		ProductElement<ComplexNumber> d11 = ring.subtract(g.entry(2, 2), ring.multiply(l10, adjoint(layer, l10), d00));
		return new Pair<>(twoByTwo(ring.one(), ring.zero(), l10, ring.one()),
				twoByTwo(d00, ring.zero(), ring.zero(), d11));
	}

	private FalconTree ffldl(int layer, Matrix<ProductElement<ComplexNumber>> g) {
		Pair<Matrix<ProductElement<ComplexNumber>>, Matrix<ProductElement<ComplexNumber>>> ldl = ldl(layer, g);
		FalconTree node = new FalconTree();
		node.innerValue = ldl.getFirst().entry(2, 1);
		if (layer == 1) {
			node.left = new FalconTree();
			node.left.leafValue = r.multiply(sigma, r.inverseSqrt(ldl.getSecond().entry(1, 1).get(0).realPart()));
			node.right = new FalconTree();
			node.right.leafValue = r.multiply(sigma, r.inverseSqrt(ldl.getSecond().entry(2, 2).get(0).realPart()));
			return node;
		}
		Vector<ProductElement<ComplexNumber>> d0 = splitFFT(layer, ldl.getSecond().entry(1, 1));
		Vector<ProductElement<ComplexNumber>> d1 = splitFFT(layer, ldl.getSecond().entry(2, 2));
		Matrix<ProductElement<ComplexNumber>> g0 = twoByTwo(d0.get(1), d0.get(2), adjoint(layer - 1, d0.get(2)),
				d0.get(1));
		Matrix<ProductElement<ComplexNumber>> g1 = twoByTwo(d1.get(1), d1.get(2), adjoint(layer - 1, d1.get(2)),
				d1.get(1));
		node.left = ffldl(layer - 1, g0);
		node.right = ffldl(layer - 1, g1);
		return node;
	}

	private NTT<PFE> reduceAndTransform(NFE t) {
		PrimeField fp = PrimeField.getPrimeField(prime);
		UnivariatePolynomialRing<PFE> ring = fp.getUnivariatePolynomialRing();
		return ntt.fromPolynomial(
				ring.getEmbedding(t.asPolynomial(), new FunctionMathMap<>((Fraction s) -> fp.getFraction(s))));
	}

	private NFE lift(NTT<PFE> t) {
		Rationals q = Rationals.q();
		PrimeField fp = PrimeField.getPrimeField(prime);
		return nf.fromPolynomial(q.getUnivariatePolynomialRing().getEmbedding(ntt.asPolynomial(t),
				new FunctionMathMap<PFE, Fraction>((PFE s) -> q.getInteger(fp.liftToInteger(s)))));
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
		return nf.fromPolynomial(q.getUnivariatePolynomialRing().getPolynomial(coeffs));
	}

	@Override
	public FalconPublicKey createPublicKey(FalconPrivateKey privateKey) {
		return privateKey.publicKey;
	}

	@Override
	public FalconPrivateKey createPrivateKey() {
		FalconPrivateKey key = ntruKeyGen();
		Matrix<ProductElement<ComplexNumber>> b = key.b;
		FreeModule<ProductElement<ComplexNumber>> module = new FreeModule<>(fftTower.get(log2degree), 2);
		MatrixAlgebra<ProductElement<ComplexNumber>> algebra = module.matrixAlgebra();
		Matrix<ProductElement<ComplexNumber>> fftBt = algebra.transpose(Matrix
				.mapMatrix(new FunctionMathMap<>((ProductElement<ComplexNumber> s) -> adjoint(log2degree, s)), b));
		Matrix<ProductElement<ComplexNumber>> g = algebra.multiply(b, fftBt);
		key.root = ffldl(log2degree, g);
		key.publicKey = new FalconPublicKey(
				lift(ntt.divideChecked(reduceAndTransform(key.g), reduceAndTransform(key.f))));
		return null;
	}

	@Override
	public FalconSignature sign(byte[] message, FalconPrivateKey privateKey) {
		byte[] random = new byte[40];
		this.random.nextBytes(random);
		ProductRing<ComplexNumber, Complex> ring = fftTower.get(log2degree);
		ProductElement<ComplexNumber> c = asFFT(hashToPoint(ByteArray.concatenate(random, message)));
		ProductElement<ComplexNumber> qInv = ring.inverse(ring.getInteger(prime));
		FreeModule<ProductElement<ComplexNumber>> module = new FreeModule<>(ring, 2);
		Vector<ProductElement<ComplexNumber>> t = new Vector<>(ring.multiply(-1, qInv, c, asFFT(privateKey.capitalF)),
				ring.multiply(qInv, c, asFFT(privateKey.f)));
		Vector<ProductElement<ComplexNumber>> s;
		do {
			Vector<ProductElement<ComplexNumber>> intSample = asComplexVector(
					fastFourierSampler(log2degree, t.get(1), t.get(2), privateKey.root));
			s = module.matrixAlgebra().multiply(module.matrixAlgebra().transpose(privateKey.b),
					module.subtract(t, intSample));
		} while (norm(s).compareTo(limit) > 0);
		NFE signature = fromFFT(fft.fromProductElement(s.get(2)));
		return new FalconSignature(random, signature);
	}

	@Override
	public boolean verify(byte[] message, FalconSignature signature, FalconPublicKey publicKey) {
		NFE c = hashToPoint(ByteArray.concatenate(signature.random.array(), message));
		NFE s1 = lift(ntt.subtract(reduceAndTransform(c),
				ntt.multiply(reduceAndTransform(publicKey.h), reduceAndTransform(signature.signature))));
		return innerProduct(log2degree, s1, signature.signature).compareTo(limit) <= 0;
	}
}
