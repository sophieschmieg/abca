package cryptography;

import static org.junit.jupiter.api.Assertions.assertTrue;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.junit.jupiter.api.Test;

import cryptography.Falcon.FalconPrivateKey;
import cryptography.Falcon.FalconPublicKey;
import cryptography.Falcon.FalconSignature;
import cryptography.FalconTest.AttackLatticeWrapper.AttackLattice;
import fields.finitefields.PrimeField;
import fields.finitefields.PrimeField.PFE;
import fields.floatingpoint.Complex;
import fields.floatingpoint.Complex.ComplexNumber;
import fields.floatingpoint.FiniteRealVectorSpace;
import fields.floatingpoint.Reals;
import fields.floatingpoint.Reals.Real;
import fields.helper.AbstractLattice;
import fields.helper.NTTRing;
import fields.helper.NTTRing.NTT;
import fields.integers.Integers;
import fields.integers.Integers.IntE;
import fields.integers.Rationals;
import fields.integers.Rationals.Fraction;
import fields.interfaces.MathMap;
import fields.numberfields.NumberField;
import fields.numberfields.NumberField.NFE;
import fields.vectors.FreeModule;
import fields.vectors.Matrix;
import fields.vectors.Vector;
import util.FunctionMathMap;
import util.Pair;

class FalconTest {

	@Test
	void test() throws IOException {
		Falcon falcon = new Falcon(9);
		FalconPrivateKey privateKey = falcon.createPrivateKey();
		FalconPublicKey publicKey = falcon.createPublicKey(privateKey);
		FalconSignature signature = falcon.sign("I love my wife".getBytes(), privateKey);
		assertTrue(falcon.verify("I love my wife".getBytes(), signature, publicKey));
	}

	@Test
	void testInsecureSign() {
		Integers z = Integers.z();
		Rationals q = Rationals.q();
		Falcon falcon = new Falcon(9);
		FalconPrivateKey privateKey = falcon.createPrivateKey();
		FalconPublicKey publicKey = falcon.createPublicKey(privateKey);
		Pair<ByteArray, NTT<ComplexNumber>> hashed = falcon.hash("I love my wife".getBytes());
		NumberField nf = falcon.getNumberField(9);
		NFE c = falcon.fromFFT(hashed.getSecond());
		List<IntE> initialList = new ArrayList<>();
		initialList.addAll(
				Vector.mapVector(new FunctionMathMap<>((Fraction t) -> t.asInteger()), nf.asVector(c)).asList());
		initialList.addAll(Vector
				.mapVector(new FunctionMathMap<>((Fraction t) -> t.asInteger()), nf.asVector(nf.zero())).asList());
		Vector<IntE> initial = new Vector<>(initialList);
		MathMap<Fraction, IntE> reductionMap = new FunctionMathMap<>((Fraction t) -> t.asInteger());
		Matrix<IntE> gMatrix = Matrix.mapMatrix(reductionMap, nf.asMatrix(nf.negative(privateKey.getG())));
		Matrix<IntE> fMatrix = Matrix.mapMatrix(reductionMap, nf.asMatrix(privateKey.getF()));
		Matrix<IntE> capitalGMatrix = Matrix.mapMatrix(reductionMap,
				nf.asMatrix(nf.negative(privateKey.getCapitalG())));
		Matrix<IntE> capitalFMatrix = Matrix.mapMatrix(reductionMap, nf.asMatrix(privateKey.getCapitalF()));
		List<List<Matrix<IntE>>> baseMatrix = new ArrayList<>();
		baseMatrix.add(new ArrayList<>());
		baseMatrix.add(new ArrayList<>());
		baseMatrix.get(0).add(fMatrix);
		baseMatrix.get(0).add(capitalFMatrix);
		baseMatrix.get(1).add(gMatrix);
		baseMatrix.get(1).add(capitalGMatrix);
		Reals r = Reals.r(64);
		MathMap<IntE, Real> embedding = new FunctionMathMap<>((IntE t) -> r.getInteger(t));
		MathMap<Real, IntE> round = new FunctionMathMap<>((Real t) -> t.round());
		Matrix<IntE> baseSystem = Matrix.fromBlockMatrix(baseMatrix);
		Matrix<Real> realBaseSystem = Matrix.mapMatrix(embedding, baseSystem);
		Vector<Real> inBase = realBaseSystem.getModule(r).solve(realBaseSystem, Vector.mapVector(embedding, initial));
		Vector<IntE> closest = baseSystem.getModule(z).multiply(baseSystem, Vector.mapVector(round, inBase));
		Vector<Fraction> solution = Vector.mapVector(q.getEmbeddingMap(),
				new Vector<>(closest.asList().subList(nf.degree(), 2 * nf.degree())));
		NFE s1 = nf.fromVector(solution);
		FalconSignature signature = falcon.new FalconSignature(hashed.getFirst().array(), s1);
		assertTrue(falcon.verify("I love my wife".getBytes(), signature, publicKey));
		// s0 + h*s1 = c
		// s0*g + f*s1 = c*g
	}

	private Real size(Falcon falcon, int log2degree, Vector<NFE> t) {
		Reals r = Reals.r(64);
		NumberField nf = falcon.getNumberField(log2degree);
		NFE f = t.get(1);
		NFE g = t.get(2);
		NTTRing<ComplexNumber, Complex> ring = falcon.getFFTRing(log2degree);
		NTT<ComplexNumber> denominator = ring.add(
				ring.multiply(falcon.asFFT(f), falcon.asFFT(falcon.adjoint(log2degree, f))),
				ring.multiply(falcon.asFFT(g), falcon.asFFT(falcon.adjoint(log2degree, g))));
		Real gamma = r.max(falcon.normSqr(log2degree, new Vector<>(f, nf.negative(g))),
				falcon.normSqr(new Vector<>(
						ring.divide(falcon.asFFT(nf.multiply(falcon.getPrime(), falcon.adjoint(log2degree, f))),
								denominator),
						ring.divide(falcon.asFFT(nf.multiply(falcon.getPrime(), falcon.adjoint(log2degree, g))),
								denominator))));
		return gamma;
	}

	class AttackLatticeWrapper {
		private Falcon falcon;
		private FiniteRealVectorSpace space;
		private FreeModule<NFE> module;

		public AttackLatticeWrapper(Falcon falcon) {
			this.falcon = falcon;
			this.space = new FiniteRealVectorSpace(Reals.r(64), 2 * falcon.getDegree());
			this.module = new FreeModule<>(falcon.getNumberField(falcon.getLog2Degree()), 2);
		}

		public AttackLattice newLattice(List<Vector<NFE>> generators) {
			return new AttackLattice(generators);
		}

		class AttackLattice extends AbstractLattice<Vector<NFE>, Real, Vector<Real>> {

			public AttackLattice(List<Vector<NFE>> generators) {
				super(space, generators, false, 0, 0.75, true);
			}

			@Override
			public Vector<Real> embedding(Vector<NFE> t) {
				List<ComplexNumber> fft = new ArrayList<>();
				fft.addAll(falcon.getFFTRing(falcon.getLog2Degree()).asList(falcon.asFFT(t.get(1))).subList(0,
						falcon.getDegree() / 2));
				fft.addAll(falcon.getFFTRing(falcon.getLog2Degree()).asList(falcon.asFFT(t.get(2))).subList(0,
						falcon.getDegree() / 2));
				List<Real> result = new ArrayList<>();
				for (ComplexNumber s : fft) {
					result.add(s.realPart());
					result.add(s.complexPart());
				}
				return new Vector<>(result);
			}

			@Override
			public Vector<NFE> zero() {
				return module.zero();
			}

			@Override
			public Vector<NFE> add(Vector<NFE> s1, Vector<NFE> s2) {
				return module.add(s1, s2);
			}

			@Override
			public Vector<NFE> negative(Vector<NFE> s) {
				return module.negative(s);
			}

			@Override
			public Vector<NFE> scalarMultiply(IntE t, Vector<NFE> s) {
				return module.scalarMultiply(falcon.getNumberField(falcon.getLog2Degree()).getInteger(t), s);
			}

			@Override
			public Exactness exactness() {
				return Exactness.EXACT;
			}

		}
	}

	@Test
	void testPrehashedKeyRecovery() {
		int log2degree = 6;
		Falcon falcon = new Falcon(log2degree);
		FalconPrivateKey privateKey = falcon.createPrivateKey();
		FalconPublicKey publicKey = falcon.createPublicKey(privateKey);
		NTTRing<ComplexNumber, Complex> fft = falcon.getFFTRing(log2degree);
		NTTRing<PFE, PrimeField> ntt = falcon.getNTTRing();
		Reals r = Reals.r(64);
		Real limit = r.multiply(1.17 * 1.17, r.getInteger(falcon.getPrime()));
		List<Vector<NFE>> basis = new ArrayList<>();
		AttackLatticeWrapper wrapper = new AttackLatticeWrapper(falcon);
		while (true) {
			FalconSignature signature = falcon.prehashSign(privateKey, new byte[40], fft.zero());
			NFE upper = falcon.lift(ntt.negative(ntt.multiply(falcon.reduceAndTransform(publicKey.getH()),
					falcon.reduceAndTransform(signature.getSignature()))));
			NFE lower = signature.getSignature();
			Vector<NFE> current = new Vector<>(upper, lower);
			System.out.println("Basis Size: " + basis.size());
			basis.add(0, current);
			if (basis.size() <= 2 * falcon.getDegree()) {
				continue;
			}
			basis = basis.subList(0, 2 * falcon.getDegree());
			AttackLattice lattice = wrapper.newLattice(basis);
			basis.clear();
			basis.addAll(lattice.getModuleGenerators());
			Vector<NFE> smallest = basis.get(0);
			NFE g = smallest.get(1);
			NFE f = smallest.get(2);
			NTT<ComplexNumber> denominator = fft.add(
					fft.multiply(falcon.asFFT(f), falcon.asFFT(falcon.adjoint(log2degree, f))),
					fft.multiply(falcon.asFFT(g), falcon.asFFT(falcon.adjoint(log2degree, g))));
			Real gamma = r
					.max(falcon.normSqr(log2degree, new Vector<>(f, falcon.getNumberField(log2degree).negative(g))),
							falcon.normSqr(
									new Vector<>(
											fft.divide(
													falcon.asFFT(falcon.getNumberField(log2degree).multiply(
															falcon.getPrime(), falcon.adjoint(log2degree, f))),
													denominator),
											fft.divide(
													falcon.asFFT(falcon.getNumberField(log2degree).multiply(
															falcon.getPrime(), falcon.adjoint(log2degree, g))),
													denominator))));
			System.out.println(gamma);
			if (gamma.compareTo(limit) < 0) {
				break;
			}
		}
		Vector<NFE> smallest = basis.get(0);
		NFE g = smallest.get(1);
		NFE f = falcon.getNumberField(log2degree).negative(smallest.get(2));
		FalconPrivateKey sk = falcon.createPrivateKeyFromGiven(f, g);
		System.out.println("Original Public Key: " + publicKey.getH());
		System.out.println("Forgery Public Key: " + falcon.createPublicKey(sk).getH());
		byte[] msg = "forge this!".getBytes();
		FalconSignature sig = falcon.sign(msg, sk);
		if (falcon.verify(msg, sig, publicKey)) {
			System.out.println("Forgery successful!");
		} else {
			System.out.println("Forgery not successful!");		
		}
	}
}
