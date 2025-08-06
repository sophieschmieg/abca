package cryptography;

import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.security.SecureRandom;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import org.junit.jupiter.api.Test;

import cryptography.Dilithium.DilithiumPrivateKey;
import cryptography.Dilithium.DilithiumPublicKey;
import cryptography.Dilithium.DilithiumSignature;
import fields.finitefields.PrimeField;
import fields.finitefields.PrimeField.PFE;
import fields.helper.NTTRing;
import fields.helper.NTTRing.NTT;
import fields.integers.Integers;
import fields.interfaces.UnivariatePolynomial;
import fields.interfaces.UnivariatePolynomialRing;
import fields.vectors.FreeModule;
import fields.vectors.Matrix;
import fields.vectors.MatrixModule;
import fields.vectors.Vector;

class DilithiumTest {

	@Test
	void test() {
		Dilithium d3 = Dilithium.DILITHIUM3;
		DilithiumPrivateKey privateKey = d3.createPrivateKey();
		DilithiumPublicKey publicKey = d3.createPublicKey(privateKey);
		byte[] message = new byte[100];
		new SecureRandom().nextBytes(message);
		DilithiumSignature signature = d3.sign(message, privateKey);
		assertTrue(d3.verify(message, signature, publicKey));
		assertFalse(d3.verify(new byte[100], signature, publicKey));
	}

	private class XorKeyInfo {
		private int[][] xorS;
		private int[][] xorE;
		private DilithiumPublicKey publicKey1;
		private Vector<NTT<PFE>> uncompressedPublicKey1;
		private DilithiumPublicKey publicKey2;
		private Vector<NTT<PFE>> uncompressedPublicKey2;

		public XorKeyInfo(Dilithium dilithium, DilithiumPrivateKey key1, DilithiumPrivateKey key2) {
			Integers z = Integers.z();
			NTTRing<PFE, PrimeField> ntt = dilithium.getNTTRing();
			this.publicKey1 = dilithium.createPublicKey(key1);
			this.uncompressedPublicKey1 = dilithium.getFreeModule().add(
					dilithium.getMatrixModule().multiply(dilithium.expandMatrix(publicKey1.getRho()), key1.getS()),
					key1.getE());
			this.publicKey2 = dilithium.createPublicKey(key2);
			this.uncompressedPublicKey2 = dilithium.getFreeModule().add(
					dilithium.getMatrixModule().multiply(dilithium.expandMatrix(publicKey2.getRho()), key2.getS()),
					key2.getE());
			int eta = dilithium.getEta();
			this.xorS = new int[dilithium.getCols()][dilithium.getDegree()];
			for (int i = 0; i < dilithium.getCols(); i++) {
				UnivariatePolynomial<PFE> coeff1 = ntt.asPolynomial(key1.getS().get(i + 1));
				UnivariatePolynomial<PFE> coeff2 = ntt.asPolynomial(key2.getS().get(i + 1));
				for (int j = 0; j < dilithium.getDegree(); j++) {
					this.xorS[i][j] = (eta
							- z.centeredLift(coeff1.univariateCoefficient(j), Dilithium.PRIME).intValueExact())
							^ (eta - z.centeredLift(coeff2.univariateCoefficient(j), Dilithium.PRIME).intValueExact());
				}
			}
			this.xorE = new int[dilithium.getRows()][dilithium.getDegree()];
			for (int i = 0; i < dilithium.getRows(); i++) {
				UnivariatePolynomial<PFE> coeff1 = ntt.asPolynomial(key1.getE().get(i + 1));
				UnivariatePolynomial<PFE> coeff2 = ntt.asPolynomial(key2.getE().get(i + 1));
				for (int j = 0; j < dilithium.getDegree(); j++) {
					this.xorE[i][j] = (eta
							- z.centeredLift(coeff1.univariateCoefficient(j), Dilithium.PRIME).intValueExact())
							^ (eta - z.centeredLift(coeff2.univariateCoefficient(j), Dilithium.PRIME).intValueExact());
				}
			}
		}
	}

	private Matrix<PFE> toIntMatrix(Dilithium dilithium, Matrix<NTT<PFE>> matrix) {
		PrimeField fp = PrimeField.getPrimeField(Dilithium.PRIME);
		UnivariatePolynomialRing<PFE> polynomials = fp.getUnivariatePolynomialRing();
		NTTRing<PFE, PrimeField> ntt = dilithium.getNTTRing();
		MatrixModule<NTT<PFE>> mm = matrix.getModule(ntt);
		FreeModule<NTT<PFE>> free = mm.domain();
		List<Vector<PFE>> result = new ArrayList<>();
		for (int i = 0; i < matrix.columns(); i++) {
			for (int j = 0; j < dilithium.getDegree(); j++) {
				Vector<NTT<PFE>> basisVector = free.scalarMultiply(ntt.fromPolynomial(polynomials.getVarPower(j)),
						free.getUnitVector(i + 1));
				Vector<NTT<PFE>> mapped = mm.multiply(matrix, basisVector);
				List<PFE> resultVector = new ArrayList<>();
				for (int k = 0; k < matrix.rows(); k++) {
					UnivariatePolynomial<PFE> mappedPolynomial = ntt.asPolynomial(mapped.get(k + 1));
					for (int l = 0; l < dilithium.getDegree(); l++) {
						resultVector.add(mappedPolynomial.univariateCoefficient(l));
					}
				}
				result.add(new Vector<>(resultVector));
			}
		}
		return Matrix.fromColumns(result);
	}

	private Vector<PFE> toIntVector(Dilithium dilithium, Vector<NTT<PFE>> t) {
		List<PFE> result = new ArrayList<>();
		for (NTT<PFE> coeff : t.asList()) {
			UnivariatePolynomial<PFE> poly = dilithium.getNTTRing().asPolynomial(coeff);
			for (int i = 0; i < dilithium.getDegree(); i++) {
				result.add(poly.univariateCoefficient(i));
			}
		}
		return new Vector<>(result);
	}

	private DilithiumPrivateKey xorKeyRecovery(Dilithium dilithium, List<XorKeyInfo> xors) {
		PrimeField fp = PrimeField.getPrimeField(Dilithium.PRIME);
		int bits = (int) Math.ceil(Math.log(2 * dilithium.getEta() + 1) / Math.log(2));
		List<List<Matrix<PFE>>> adjustedBlocks = new ArrayList<>();
		List<List<Matrix<PFE>>> blocks = new ArrayList<>();
		List<PFE> rhs = new ArrayList<>();
		blocks.add(new ArrayList<>());
		blocks.get(0).add(toIntMatrix(dilithium, dilithium.expandMatrix(xors.get(0).publicKey1.getRho())));
		blocks.get(0).add(toIntMatrix(dilithium, dilithium.getFreeModule().matrixAlgebra().one()));
		rhs.addAll(toIntVector(dilithium, xors.get(0).uncompressedPublicKey1).asList());
		List<Matrix<PFE>> diagonalBlocks = new ArrayList<>();
		MatrixModule<PFE> mm = new MatrixModule<>(fp, dilithium.getRows() * dilithium.getDegree() * (xors.size() + 1),
				2 * dilithium.getCols() * dilithium.getDegree() * bits);
		List<PFE> constant = new ArrayList<>();
		List<PFE> rhsAdjustment = new ArrayList<>();
		for (int i = 0; i < dilithium.getRows() + dilithium.getCols(); i++) {
			for (int j = 0; j < dilithium.getDegree(); j++) {
				List<PFE> row = new ArrayList<>();
				constant.add(fp.getElement(dilithium.getEta()));
				for (int k = 0; k < bits; k++) {
					row.add(fp.getElement(-(1 << k)));
				}
				diagonalBlocks.add(Matrix.fromRows(Collections.singletonList(new Vector<>(row))));
			}
		}
		Matrix<PFE> blockMatrix = Matrix.fromBlockMatrix(blocks);
		Matrix<PFE> adjustment = Matrix.fromBlockDiagonalMatrix(diagonalBlocks, fp.zero());
		adjustedBlocks.add(Collections.singletonList(mm.multiply(blockMatrix, adjustment)));
		rhsAdjustment.addAll(mm.multiply(blockMatrix, new Vector<>(constant)).asList());
		for (XorKeyInfo info : xors) {
			blocks.clear();
			diagonalBlocks.clear();
			constant.clear();
			blocks.add(new ArrayList<>());
			blocks.get(0).add(toIntMatrix(dilithium, dilithium.expandMatrix(info.publicKey2.getRho())));
			blocks.get(0).add(toIntMatrix(dilithium, dilithium.getFreeModule().matrixAlgebra().one()));
			rhs.addAll(toIntVector(dilithium, info.uncompressedPublicKey2).asList());
			for (int i = 0; i < dilithium.getRows() + dilithium.getCols(); i++) {
				for (int j = 0; j < dilithium.getDegree(); j++) {
					int c = dilithium.getEta();
					List<PFE> row = new ArrayList<>();
					for (int k = 0; k < bits; k++) {
						boolean xor = ((1 << k) & (i < dilithium.getCols() ? info.xorS[i][j]
								: info.xorE[i - dilithium.getRows()][j])) != 0;
						row.add(fp.getElement((xor ? 1 : -1) * (1 << k)));
						if (xor) {
							c -= 1 << k;
						}
					}
					constant.add(fp.getElement(c));
					diagonalBlocks.add(Matrix.fromRows(Collections.singletonList(new Vector<>(row))));
				}
			}
			blockMatrix = Matrix.fromBlockMatrix(blocks);
			adjustment = Matrix.fromBlockDiagonalMatrix(diagonalBlocks, fp.zero());
			adjustedBlocks.add(Collections.singletonList(mm.multiply(blockMatrix, adjustment)));
			rhsAdjustment.addAll(mm.multiply(blockMatrix, new Vector<>(constant)).asList());
		}
		Matrix<PFE> full = Matrix.fromBlockMatrix(adjustedBlocks);
		Vector<PFE> rhsVector = new Vector<>(rhs);
		rhsVector = mm.codomain().subtract(rhsVector, new Vector<>(rhsAdjustment));
		Vector<PFE> solved = mm.solve(full, rhsVector);
		System.out.println(solved);
		List<NTT<PFE>> s = new ArrayList<>();
		List<NTT<PFE>> e = new ArrayList<>();
		for (int i = 0; i < dilithium.getCols(); i++) {
			List<PFE> coeffs = new ArrayList<>();
			for (int j = 0; j < dilithium.getDegree(); j++) {
				int coeff = 0;
				for (int k = 0; k < bits; k++) {
					coeff += (1 << k) * solved.get(i * dilithium.getDegree() * bits + j * bits + k + 1).getValue()
							.intValueExact();
				}
				coeffs.add(fp.getInteger(dilithium.getEta() - coeff));
			}
			s.add(dilithium.getNTTRing().fromPolynomial(fp.getUnivariatePolynomialRing().getPolynomial(coeffs)));
		}
		for (int i = 0; i < dilithium.getRows(); i++) {
			List<PFE> coeffs = new ArrayList<>();
			for (int j = 0; j < dilithium.getDegree(); j++) {
				int coeff = 0;
				for (int k = 0; k < bits; k++) {
					coeff += (1 << k) * solved.get(dilithium.getCols() * dilithium.getDegree() * bits
							+ i * dilithium.getDegree() * bits + j * bits + k + 1).getValue().intValueExact();
				}
				coeffs.add(fp.getInteger(dilithium.getEta() - coeff));
			}
			e.add(dilithium.getNTTRing().fromPolynomial(fp.getUnivariatePolynomialRing().getPolynomial(coeffs)));
		}
		return dilithium.createPrivateKeyFromParameters(xors.get(0).publicKey1.getRho(), new byte[32], new Vector<>(s),
				new Vector<>(e));
	}

	@Test
	void testXorKeyRecovery() {
		// 29.615.817 s
		PrimeField fp = PrimeField.getPrimeField(Dilithium.PRIME);
		Dilithium dilithium = Dilithium.DILITHIUM2; /*
													 * new Dilithium(Dilithium.PRIME, 8, fp.power(fp.getElement(1753),
													 * 256 / 8).getValue().intValueExact(), 4, 4, 2, 13, 39, 17,
													 * (Dilithium.PRIME - 1) / 88, 80, 4, Sha3.SHAKE_128,
													 * Sha3.SHAKE_256, Sha3.SHAKE_256);
													 */
		List<XorKeyInfo> xors = new ArrayList<>();
		DilithiumPrivateKey key = dilithium.createPrivateKey();
		int bits = (int) Math.ceil(Math.log(2 * dilithium.getEta() + 1) / Math.log(2));
		for (int i = 0; i < 2 * bits; i++) {
			DilithiumPrivateKey key2 = dilithium.createPrivateKey();
			xors.add(new XorKeyInfo(dilithium, key, key2));
		}
		DilithiumPrivateKey recovered = xorKeyRecovery(dilithium, xors);
		System.out.println(key.getS());
		System.out.println(recovered.getS());
		System.out.println(key.getE());
		System.out.println(recovered.getE());
		// DilithiumSignature signature = dilithium.sign("This is a test
		// message".getBytes(), recovered);
		// assertTrue(dilithium.verify("This is a test message".getBytes(), signature,
		// dilithium.createPublicKey(key)));
	}
}
