package fields.helper;

import java.util.ArrayList;
import java.util.List;

import fields.interfaces.AlgebraicExtensionElement;
import fields.interfaces.Element;
import fields.interfaces.Field.Extension;
import fields.interfaces.FieldExtension;
import fields.interfaces.MathMap;
import fields.interfaces.MathSet.Exactness;
import fields.interfaces.UnivariatePolynomial;
import fields.interfaces.UnivariatePolynomialRing;
import fields.vectors.FiniteVectorSpace;
import fields.vectors.Matrix;
import fields.vectors.MatrixAlgebra;
import fields.vectors.SubVectorSpace;
import fields.vectors.Vector;
import util.Pair;

public class FieldEmbedding<T extends Element<T>, S extends AlgebraicExtensionElement<T, S>, Ext extends FieldExtension<T, S, Ext>> {
	private Ext field;
	private Ext embeddedField;
	private S embeddedAlpha;
	private S generator;
	private MathMap<S, S> embedding;
	private MathMap<S, Vector<S>> asVectorMap;
	private UnivariatePolynomial<S> minimalPolynomial;
	private Matrix<T> alphaToAlphaBetaBaseChange;
	private FiniteVectorSpace<S> asVectorSpace;
	private MatrixAlgebra<S> algebra;
	private SubVectorSpace<T, S> asSubVectorSpace;

	@SuppressWarnings("unchecked")
	public FieldEmbedding(Ext field, Ext embeddedField, MathMap<S, S> embedding) {
		this.field = field;
		this.embeddedField = embeddedField;
		if (!this.field.getBaseField().equals(this.embeddedField.getBaseField())) {
			throw new ArithmeticException("Different base fields!");
		}
		this.embedding = embedding;
		this.embeddedAlpha = embedding.evaluate(embeddedField.alpha());
		UnivariatePolynomialRing<S> polynomials = this.field.getUnivariatePolynomialRing();
		if (!polynomials.evaluate(polynomials.getEmbedding(embeddedField.minimalPolynomial(), field.getEmbeddingMap()),
				embeddedAlpha).equals(field.zero())) {
			throw new ArithmeticException("Not properly embedded!");
		}
	}

	public FieldEmbedding(Ext field, Ext embeddedField, S embeddedAlpha) {
		this(field, embeddedField, embeddedAlpha, field.alpha());
	}

	@SuppressWarnings("unchecked")
	public FieldEmbedding(Ext field, Ext embeddedField, S embeddedAlpha, S generator) {
		this.field = field;
		this.embeddedField = embeddedField;
		if (!this.field.getBaseField().equals(this.embeddedField.getBaseField())) {
			throw new ArithmeticException("Different base fields!");
		}
		this.embeddedAlpha = embeddedAlpha;
		this.generator = generator;
		this.embedding = new MathMap<S, S>() {
			private UnivariatePolynomialRing<S> polynomials = FieldEmbedding.this.field.getUnivariatePolynomialRing();

			@Override
			public S evaluate(S t) {
				return polynomials.evaluate(
						polynomials.getEmbedding(t.asPolynomial(), FieldEmbedding.this.field.getEmbeddingMap()),
						FieldEmbedding.this.embeddedAlpha);
			}

		};
		this.asVectorMap = new MathMap<S, Vector<S>>() {
			@Override
			public Vector<S> evaluate(S t) {
				return asVector(t);
			}
		};
		UnivariatePolynomialRing<S> polynomials = this.field.getUnivariatePolynomialRing();
		if (field.exactness() == Exactness.EXACT && !polynomials.evaluate(polynomials.getEmbedding(embeddedField.minimalPolynomial(), field.getEmbeddingMap()),
				embeddedAlpha).equals(field.zero())) {
			throw new ArithmeticException("Not properly embedded!");
		}
	}

	public FieldEmbedding(Extension<S, T, S, Ext> extension, Ext base) {
		this(extension.extension(), base, extension.embeddingMap().evaluate(base.alpha()));
	}

	public FieldEmbedding(Ext field, Ext embeddedField) {
		this(field, embeddedField,
				field.roots(field.getUnivariatePolynomialRing().getEmbedding(embeddedField.minimalPolynomial(),
						embeddedField.getEmbeddingMap())).keySet().iterator().next());
	}

	public FieldEmbedding(Ext field) {
		this(field, field, field.alpha());
	}

	public FieldEmbedding(FieldEmbedding<T, S, Ext> first, FieldEmbedding<T, S, Ext> second) {
		this(second.getField(), first.getEmbeddedField(), second.getEmbedding(first.embeddedAlpha));
		if (!first.getField().equals(second.getEmbeddedField())) {
			throw new ArithmeticException("Not a field extension tower");
		}
	}
	
	@Override
	public String toString() {
		return field + " / " + embeddedField;
	}

	public Ext getField() {
		return field;
	}

	public Ext getEmbeddedField() {
		return embeddedField;
	}

	public S getEmbedding(S s) {
		return embedding.evaluate(s);
	}

	public MathMap<S, S> getEmbeddingMap() {
		return embedding;
	}

	public MathMap<S, Vector<S>> asVectorMap() {
		return asVectorMap;
	}
	
	public S getEmbeddedAlpha() {
		return embeddedAlpha;
	}
	
	public S getGenerator() {
		return generator;
	}

	public int relativeDegree() {
		return field.degree() / embeddedField.degree();
	}

	public Vector<S> asVector(S s) {
		List<T> alphaBetaVector = asAlphaBetaVector(s).asList();
		List<S> asList = new ArrayList<>();
		for (int i = 0; i < relativeDegree(); i++) {
			asList.add(embeddedField.fromVector(new Vector<>(
					alphaBetaVector.subList(i * embeddedField.degree(), (i + 1) * embeddedField.degree()))));
		}
		return new Vector<>(asList);
	}
	
	public S fromVector(Vector<S> s) {
		return fromPolynomial(embeddedField.getUnivariatePolynomialRing().getPolynomial(s.asList()));
	}
	
	public UnivariatePolynomial<S> asPolynomial(S s) {
		return embeddedField.getUnivariatePolynomialRing().getPolynomial(asVector(s).asList());
	}

	@SuppressWarnings("unchecked")
	public S fromPolynomial(UnivariatePolynomial<S> polynomial) {
		UnivariatePolynomialRing<S> ring = field.getUnivariatePolynomialRing();
		return ring.evaluate(ring.getEmbedding(polynomial, embedding), generator);
	}

	public Matrix<T> alphaToAlphaBetaBaseChange() {
		if (alphaToAlphaBetaBaseChange == null) {
			List<Vector<T>> basis = new ArrayList<>();
			for (int i = 0; i < relativeDegree(); i++) {
				for (int j = 0; j < embeddedField.degree(); j++) {
					S generator = field.multiply(field.power(embeddedAlpha, j), field.power(this.generator, i));
					basis.add(field.asVector(generator));
				}
			}
			Matrix<T> matrix = Matrix.fromColumns(basis);
			this.alphaToAlphaBetaBaseChange = this.field.matrixAlgebra().inverse(matrix);
		}
		return alphaToAlphaBetaBaseChange;
	}

	public Vector<T> asAlphaBetaVector(S s) {
		return field.matrixAlgebra().multiply(alphaToAlphaBetaBaseChange(), field.asVector(s));
	}

	public UnivariatePolynomial<S> minimalPolynomial() {
		if (minimalPolynomial == null) {
			this.minimalPolynomial = field.minimalPolynomialOver(generator, this);
		}
		return minimalPolynomial;
	}

	public UnivariatePolynomial<S> minimalPolynomial(S s) {
		return field.minimalPolynomialOver(s, this);
	}

	public Matrix<S> asMatrix(S s) {
		return field.asMatrixOver(s, this);
	}

	public FiniteVectorSpace<S> asVectorSpace() {
		if (asVectorSpace == null) {
			asVectorSpace = new FiniteVectorSpace<>(embeddedField, relativeDegree());
		}
		return asVectorSpace;
	}

	public MatrixAlgebra<S> matrixAlgebra() {
		if (algebra == null) {
			algebra = asVectorSpace().matrixAlgebra();
		}
		return algebra;
	}

	public S norm(S s) {
		return matrixAlgebra().determinant(asMatrix(s));
	}

	public S trace(S s) {
		return matrixAlgebra().trace(asMatrix(s));
	}

	public SubVectorSpace<T, S> asSubVectorSpace() {
		if (asSubVectorSpace == null) {
			List<S> generators = new ArrayList<>();
			S generator = field.one();
			for (int i = 0; i < embeddedField.degree(); i++) {
				generators.add(generator);
				generator = field.multiply(generator, embeddedAlpha);
			}
			asSubVectorSpace = new SubVectorSpace<>(field, generators);
		}
		return asSubVectorSpace;
	}

	public static class CompositionResult<T extends Element<T>, S extends AlgebraicExtensionElement<T, S>, Ext extends FieldExtension<T, S, Ext>> {
		private FieldEmbedding<T, S, Ext> composition;
		private FieldEmbedding<T, S, Ext> thisEmbedding;
		private FieldEmbedding<T, S, Ext> otherEmbedding;
		private FieldEmbedding<T, S, Ext> intersection;
		private FieldEmbedding<T, S, Ext> thisInComposition;
		private FieldEmbedding<T, S, Ext> otherInComposition;
		private FieldEmbedding<T, S, Ext> intersectionInThis;
		private FieldEmbedding<T, S, Ext> intersectionInOther;
		private FieldEmbedding<T, S, Ext> intersectionInComposition;

		protected CompositionResult(FieldEmbedding<T, S, Ext> composition, FieldEmbedding<T, S, Ext> thisEmbedding,
				FieldEmbedding<T, S, Ext> otherEmbedding, FieldEmbedding<T, S, Ext> intersection) {
			this.composition = composition;
			this.thisEmbedding = thisEmbedding;
			this.otherEmbedding = otherEmbedding;
			this.intersection = intersection;
			this.thisInComposition = new FieldEmbedding<>(composition.getEmbeddedField(),
					thisEmbedding.getEmbeddedField(), composition.asVector(thisEmbedding.embeddedAlpha).get(1));
			this.otherInComposition = new FieldEmbedding<>(composition.getEmbeddedField(),
					otherEmbedding.getEmbeddedField(), composition.asVector(otherEmbedding.embeddedAlpha).get(1));
			this.intersectionInThis = new FieldEmbedding<>(thisEmbedding.getEmbeddedField(),
					intersection.getEmbeddedField(), thisEmbedding.asVector(intersection.embeddedAlpha).get(1));
			this.intersectionInOther = new FieldEmbedding<>(otherEmbedding.getEmbeddedField(),
					intersection.getEmbeddedField(), otherEmbedding.asVector(intersection.embeddedAlpha).get(1));
			this.intersectionInComposition = new FieldEmbedding<>(composition.getEmbeddedField(),
					intersection.getEmbeddedField(), composition.asVector(intersection.embeddedAlpha).get(1));
		}

		public FieldEmbedding<T, S, Ext> getThisInCompositionEmbedding() {
			return thisInComposition;
		}

		public FieldEmbedding<T, S, Ext> getOtherInCompositionEmbedding() {
			return otherInComposition;
		}

		public FieldEmbedding<T, S, Ext> getIntersectionInThisEmbedding() {
			return intersectionInThis;
		}

		public FieldEmbedding<T, S, Ext> getIntersectionInOtherEmbedding() {
			return intersectionInOther;
		}

		public FieldEmbedding<T, S, Ext> getIntersectionInCompositionEmbedding() {
			return intersectionInComposition;
		}

		public FieldEmbedding<T, S, Ext> getComposition() {
			return composition;
		}

		public FieldEmbedding<T, S, Ext> getThisEmbedding() {
			return thisEmbedding;
		}

		public FieldEmbedding<T, S, Ext> getOtherEmbedding() {
			return otherEmbedding;
		}

		public FieldEmbedding<T, S, Ext> intersection() {
			return intersection;
		}
	}

	public CompositionResult<T, S, Ext> composeWith(FieldEmbedding<T, S, Ext> other) {
		if (!this.field.equals(other.field)) {
			throw new ArithmeticException("Not subfield of a common super field");
		}
		if (!this.embeddedField.getBaseField().equals(other.embeddedField.getBaseField())) {
			throw new ArithmeticException("No common base field");
		}
		SubVectorSpace<T, S> thisSubVectorSpace = this.asSubVectorSpace();
		SubVectorSpace<T, S> otherSubVectorSpace = other.asSubVectorSpace();
		SubVectorSpace<T, S> intersectionAsSubVectorSpace = thisSubVectorSpace.intersection(otherSubVectorSpace);
		SubVectorSpace<T, S> compositionAsSubVectorSpace = thisSubVectorSpace.add(otherSubVectorSpace);
		Pair<S, UnivariatePolynomial<T>> intersectionData = getGenerator(intersectionAsSubVectorSpace);
		S intersectionAlpha = intersectionData.getFirst();
		FieldEmbedding<T, S, Ext> intersection = new FieldEmbedding<T, S, Ext>(field,
				field.makeExtension(intersectionData.getSecond()), intersectionAlpha);
		Pair<S, UnivariatePolynomial<T>> compositionData = getGenerator(compositionAsSubVectorSpace);
		S compositionAlpha = compositionData.getFirst();
		FieldEmbedding<T, S, Ext> composition = new FieldEmbedding<T, S, Ext>(field,
				field.makeExtension(compositionData.getSecond()), compositionAlpha);
		return new CompositionResult<>(composition, this, other, intersection);
	}

	private Pair<S, UnivariatePolynomial<T>> getGenerator(SubVectorSpace<T, S> asSubVectorSpace) {
		for (S s : asSubVectorSpace.getModuleGenerators()) {
			UnivariatePolynomial<T> mipo = field.minimalPolynomial(s);
			if (mipo.degree() == asSubVectorSpace.dimension()) {
				return new Pair<>(s, mipo);
			}
		}
		while (true) {
			S s = asSubVectorSpace.getRandomElement();
			UnivariatePolynomial<T> mipo = field.minimalPolynomial(s);
			if (mipo.degree() == asSubVectorSpace.dimension()) {
				return new Pair<>(s, mipo);
			}
		}
	}

}
