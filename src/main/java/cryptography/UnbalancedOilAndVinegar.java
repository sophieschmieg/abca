package cryptography;

import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import cryptography.UnbalancedOilAndVinegar.UOVPrivateKey;
import cryptography.interfaces.ExtendedOutputFunction;
import cryptography.interfaces.SignatureScheme;
import fields.finitefields.Char2Field;
import fields.finitefields.Char2Field.C2E;
import fields.finitefields.FiniteField;
import fields.finitefields.PrimeField;
import fields.helper.AbstractElement;
import fields.interfaces.Polynomial;
import fields.interfaces.PolynomialRing;
import fields.polynomials.AbstractPolynomialRing;
import fields.polynomials.Monomial;
import fields.vectors.FiniteVectorSpace;
import fields.vectors.FreeModule;
import fields.vectors.Matrix;
import fields.vectors.MatrixModule;
import fields.vectors.Vector;
import util.MiscAlgorithms;

public class UnbalancedOilAndVinegar
		implements SignatureScheme<Vector<C2E>, Vector<Polynomial<C2E>>, UOVPrivateKey, UnbalancedOilAndVinegar> {
	public class UOVPrivateKey extends AbstractElement<UOVPrivateKey> {
		private Matrix<C2E> oilMatrix;
		private Vector<Polynomial<C2E>> vinegarForms;
		private Vector<Matrix<C2E>> baseChange;
		private Vector<Polynomial<C2E>> publicKey;

		public UOVPrivateKey(Vector<Polynomial<C2E>> publicKey, Matrix<C2E> oilMatrix,
				Vector<Polynomial<C2E>> vinegarForms, Vector<Matrix<C2E>> baseChange) {
			super();
			this.publicKey = publicKey;
			this.oilMatrix = oilMatrix;
			this.vinegarForms = vinegarForms;
			this.baseChange = baseChange;
		}

		@Override
		public int compareTo(UOVPrivateKey o) {
			int cmp = publicKey.compareTo(o.publicKey);
			if (cmp != 0) {
				return cmp;
			}
			cmp = oilMatrix.compareTo(o.oilMatrix);
			if (cmp != 0) {
				return cmp;
			}
			cmp = vinegarForms.compareTo(o.vinegarForms);
			if (cmp != 0) {
				return cmp;
			}
			return baseChange.compareTo(o.baseChange);
		}

		public Vector<Polynomial<C2E>> getPublicKey() {
			return publicKey;
		}

		public Matrix<C2E> getOilMatrix() {
			return oilMatrix;
		}

		public Vector<Polynomial<C2E>> getVinegarForms() {
			return vinegarForms;
		}

		public Vector<Matrix<C2E>> getBaseChange() {
			return baseChange;
		}

	}

	private int power;
	private int oil;
	private int vinegar;
	private PolynomialRing<C2E> vinegarPolynomialRing;
	private PolynomialRing<C2E> polynomialRing;
	private Char2Field field;
	private ExtendedOutputFunction xof;

	public static final UnbalancedOilAndVinegar uovIp = new UnbalancedOilAndVinegar(44, 112 - 44,
			Char2Field.fromLittleEndian(8, new byte[] { 0x1b }), Sha3.SHAKE_256);

	public UnbalancedOilAndVinegar(int oil, int vinegar, int power, ExtendedOutputFunction xof) {
		this(oil, vinegar, Char2Field.fromFiniteField(FiniteField.getFiniteField(PrimeField.getPrimeField(2), power)),
				xof);
	}

	public UnbalancedOilAndVinegar(int oil, int vinegar, Char2Field field, ExtendedOutputFunction xof) {
		this.power = field.degree();
		this.oil = oil;
		this.vinegar = vinegar;
		this.field = field;
		this.vinegarPolynomialRing = AbstractPolynomialRing.getPolynomialRing(field, vinegar, Monomial.GREVLEX);
		this.polynomialRing = AbstractPolynomialRing.getPolynomialRing(field, oil + vinegar, Monomial.GREVLEX);
		this.xof = xof;
	}

	public Vector<C2E> parse(byte[] msg) {
		InputStream hashed = xof.stream(msg);
		List<C2E> parsed = new ArrayList<>();
		for (int i = 0; i < oil; i++) {
			byte[] elements = new byte[MiscAlgorithms.DivRoundUp(power, 8)];
			try {
				if (hashed.read(elements) != elements.length) {
					throw new ArithmeticException("XOF ended unexpectedly");
				}
			} catch (IOException e) {
				throw new ArithmeticException("XOF ended unexpectedly");
			}
			parsed.add(field.fromBytes(elements));
		}
		return new Vector<>(parsed);
	}

	public Char2Field getField() {
		return field;
	}

	public PolynomialRing<C2E> getPolynomialRing() {
		return polynomialRing;
	}

	public Matrix<C2E> asMatrix(Polynomial<C2E> polynomial, int dim) {
		if (!polynomialRing.isHomogeneous(polynomial) || polynomial.degree() != 2) {
			throw new ArithmeticException("Not a quadratic form");
		}
		for (int i = dim; i < polynomialRing.numberOfVariables(); i++) {
			if (polynomial.degree(i + 1) != 0) {
				throw new ArithmeticException("Dimension too small!");
			}
		}
		List<List<C2E>> result = new ArrayList<>();
		for (int i = 0; i < dim; i++) {
			List<C2E> row = new ArrayList<>();
			for (int j = 0; j < i; j++) {
				row.add(field.zero());
			}
			for (int j = i; j < dim; j++) {
				int[] exp = new int[polynomialRing.numberOfVariables()];
				exp[i]++;
				exp[j]++;
				row.add(polynomial.coefficient(polynomialRing.getMonomial(exp)));
			}
			result.add(row);
		}
		return new Matrix<>(result);
	}

	public Polynomial<C2E> fromMatrix(Matrix<C2E> matrix, PolynomialRing<C2E> polynomialRing) {
		if (matrix.columns() != matrix.rows()) {
			throw new ArithmeticException("Not a square matrix!");
		}
		if (matrix.rows() != polynomialRing.numberOfVariables()) {
			throw new ArithmeticException("Wrong polynomial ring!");
		}
		Matrix<Polynomial<C2E>> embedding = Matrix.mapMatrix(polynomialRing.getEmbeddingMap(), matrix);
		MatrixModule<Polynomial<C2E>> mm = embedding.getModule(polynomialRing);
		FreeModule<Polynomial<C2E>> domain = mm.domain();
		List<Polynomial<C2E>> vars = new ArrayList<>();
		for (int i = 0; i < matrix.columns(); i++) {
			vars.add(polynomialRing.getVar(i + 1));
		}
		Vector<Polynomial<C2E>> vector = new Vector<>(vars);
		Vector<Polynomial<C2E>> eval = mm.multiply(embedding, vector);
		return domain.innerProduct(vector, eval);
	}

	@Override
	public Vector<Polynomial<C2E>> createPublicKey(UOVPrivateKey privateKey) {
		return privateKey.getPublicKey();
	}

	@Override
	public UOVPrivateKey createPrivateKey() {
		List<List<C2E>> oilSpace = new ArrayList<>();
		for (int i = 0; i < vinegar; i++) {
			List<C2E> row = new ArrayList<>();
			for (int j = 0; j < oil; j++) {
				row.add(field.getRandomElement());
			}
			oilSpace.add(row);
		}
		Matrix<C2E> oilMatrix = new Matrix<>(oilSpace);
		// B = (I O)
		// (0 I)
		// A = (A11 A12)
		// ( 0 A22)
		// B^T.A.B = (I 0).(A11 A12).(I O) = (I 0).(A11 A11.O + A12)
		// (O^T I).( 0 A22).(0 I) (O^T I).( 0 A22)
		// = (A11 A11.O + A12)
		// (O^T.A11 O^T.A11.O + O^T.A12 + A22)

		// A22 = -O^T.A11.O - O^T.A12

		// v, v^T.A11.v
		// v^T.S.x = t
		// (v^T + (O.x)^T x^T).A.(v + Ox, x) = t
		// (v + Ox)^T.A11.(v + Ox) + (v + Ox)^T.A12.x + x^T.A22.x = t
		// v^T.A11.v + x^T.O^T.A11.v + v^T.A11.O.x + x^T.O^T.A11.O.x + v^T.A12.x +
		// x^T.O^T.A12.x - x^T.O^T.A11.O.x - x^T.O^T.A12.x = t
		// v^T.A11.v + x^T.O^T.A11.v + v^T.A11.O.x + v^T.A12.x = t
		// v^T.A11.v + v^T.A11^T.O.x + v^T.A11.O.x + v^T.A12.x = t
		// v^T.((A11^T + A11).O + A12).x = t - v^T.A11.v
		List<Polynomial<C2E>> publicKey = new ArrayList<>();
		List<Polynomial<C2E>> vinegarForms = new ArrayList<>();
		List<Matrix<C2E>> baseChanges = new ArrayList<>();
		for (int k = 0; k < oil; k++) {
			List<List<C2E>> vinegarVinegar = new ArrayList<>();
			for (int i = 0; i < vinegar; i++) {
				List<C2E> row = new ArrayList<>();
				for (int j = 0; j < i; j++) {
					row.add(field.zero());
				}
				for (int j = i; j < vinegar; j++) {
					row.add(field.getRandomElement());
				}
				vinegarVinegar.add(row);
			}
			Matrix<C2E> a11 = new Matrix<>(vinegarVinegar);
			vinegarForms.add(fromMatrix(a11, vinegarPolynomialRing));
			List<List<C2E>> vinegarOil = new ArrayList<>();
			for (int i = 0; i < vinegar; i++) {
				List<C2E> row = new ArrayList<>();
				for (int j = 0; j < oil; j++) {
					row.add(field.getRandomElement());
				}
				vinegarOil.add(row);
			}
			Matrix<C2E> a12 = new Matrix<>(vinegarOil);
			MatrixModule<C2E> upperLeft = new MatrixModule<Char2Field.C2E>(field, vinegar, vinegar);
			MatrixModule<C2E> upperRight = new MatrixModule<Char2Field.C2E>(field, vinegar, oil);
			MatrixModule<C2E> lowerLeft = new MatrixModule<Char2Field.C2E>(field, oil, vinegar);
			MatrixModule<C2E> lowerRight = new MatrixModule<Char2Field.C2E>(field, oil, oil);
			Matrix<C2E> a22 = lowerRight.negative(lowerRight.multiply(upperRight.transpose(oilMatrix),
					upperRight.add(upperRight.multiply(a11, oilMatrix), a12)));
			List<List<Matrix<C2E>>> assembled = new ArrayList<>();
			assembled.add(new ArrayList<>());
			assembled.add(new ArrayList<>());
			assembled.get(0).add(a11);
			assembled.get(0).add(a12);
			assembled.get(1).add(lowerLeft.zero());
			assembled.get(1).add(a22);
			Matrix<C2E> a = Matrix.fromBlockMatrix(assembled);
			publicKey.add(fromMatrix(a, polynomialRing));
			Matrix<C2E> baseChange = upperRight
					.add(upperRight.multiply(upperLeft.add(upperLeft.transpose(a11), a11), oilMatrix), a12);
			baseChanges.add(baseChange);
		}
		return new UOVPrivateKey(new Vector<>(publicKey), oilMatrix, new Vector<>(vinegarForms),
				new Vector<>(baseChanges));

	}

	@Override
	public Vector<C2E> sign(byte[] message, UOVPrivateKey privateKey) {
		Vector<C2E> hash = parse(message);
		while (true) {
			List<C2E> saltList = new ArrayList<>();
			for (int i = 0; i < vinegar; i++) {
				saltList.add(field.getRandomElement());
			}
			Vector<C2E> salt = new Vector<>(saltList);
			Vector<C2E> eval = vinegarPolynomialRing.evaluate(privateKey.vinegarForms, salt);
			MatrixModule<C2E> oilModule = new MatrixModule<>(field, oil, oil);
			List<Vector<C2E>> rows = new ArrayList<>();
			for (int i = 0; i < oil; i++) {
				rows.add(oilModule.multiply(oilModule.transpose(privateKey.baseChange.get(i + 1)), salt));
			}
			Matrix<C2E> matrix = Matrix.fromRows(rows);
			if (oilModule.rank(matrix) != oil) {
				continue;
			}
			Vector<C2E> x = oilModule.solve(matrix, oilModule.codomain().subtract(hash, eval));
			MatrixModule<C2E> mm = privateKey.oilMatrix.getModule(field);
			List<C2E> solution = new ArrayList<>();
			solution.addAll(mm.codomain().add(salt, mm.multiply(privateKey.oilMatrix, x)).asList());
			solution.addAll(x.asList());
			return new Vector<>(solution);
		}
	}

	@Override
	public boolean verify(byte[] message, Vector<C2E> signature, Vector<Polynomial<C2E>> publicKey) {
		Vector<C2E> hash = parse(message);
		return polynomialRing.evaluate(publicKey, signature).equals(hash);
	}
}
