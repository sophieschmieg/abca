package fields.floatingpoint;

import java.util.Iterator;
import java.util.List;

import fields.floatingpoint.Complex.ComplexNumber;
import fields.floatingpoint.Reals.Real;
import fields.interfaces.ValueField;
import fields.vectors.AbstractInnerProductSpace;
import fields.vectors.FiniteVectorSpace;
import fields.vectors.MatrixAlgebra;
import fields.vectors.Vector;

public class FiniteComplexVectorSpace extends AbstractInnerProductSpace<ComplexNumber, Vector<ComplexNumber>> {
	private Complex c;
	private FiniteVectorSpace<ComplexNumber> vectorSpace;
	private int dimension;

	public FiniteComplexVectorSpace(Complex c, int dimension) {
		this.c = c;
		this.dimension = dimension;
		this.vectorSpace = new FiniteVectorSpace<>(this.c, this.dimension);
	}

	@Override
	public Exactness exactness() {
		return Exactness.FLOATING_POINT;
	}

	@Override
	public ComplexNumber innerProduct(Vector<ComplexNumber> s1, Vector<ComplexNumber> s2) {
		ComplexNumber result = c.zero();
		for (int i = 0; i < dimension; i++) {
			result = c.add(result, c.multiply(s1.get(i + 1), c.conjugate(s2.get(i + 1))));
		}
		return result;
	}

	@Override
	public ValueField<ComplexNumber> getValueField() {
		return c;
	}
	
	@Override
	public ComplexNumber asComplexNumber(ComplexNumber t) {
		return t;
	}
	
	@Override
	public ComplexNumber fromComplexNumber(ComplexNumber t) {
		return t;
	}
	
	@Override
	public ComplexNumber fromReal(Real r) {
		return c.getEmbedding(r);
	}
	
	@Override
	public ComplexNumber conjugate(ComplexNumber t) {
		return c.conjugate(t);
	}
	
	@Override
	public FiniteComplexVectorSpace withDimension(int dimension) {
		return new FiniteComplexVectorSpace(c, dimension);
	}

	@Override
	public Vector<ComplexNumber> zero() {
		return vectorSpace.zero();
	}

	@Override
	public Vector<ComplexNumber> add(Vector<ComplexNumber> s1, Vector<ComplexNumber> s2) {
		return vectorSpace.add(s1, s2);
	}

	@Override
	public Vector<ComplexNumber> negative(Vector<ComplexNumber> s) {
		return vectorSpace.negative(s);
	}

	@Override
	public Vector<ComplexNumber> scalarMultiply(ComplexNumber t, Vector<ComplexNumber> s) {
		return vectorSpace.scalarMultiply(t, s);
	}

	@Override
	public Vector<ComplexNumber> getRandomElement() {
		return vectorSpace.getRandomElement();
	}

	@Override
	public Iterator<Vector<ComplexNumber>> iterator() {
		return vectorSpace.iterator();
	}
	
	@Override
	public Vector<ComplexNumber> getUnitVector(int index) {
		return vectorSpace.getUnitVector(index);
	}

	@Override
	public List<Vector<ComplexNumber>> getBasis() {
		return vectorSpace.getBasis();
	}

	@Override
	public boolean isGeneratingModule(List<Vector<ComplexNumber>> s) {
		return vectorSpace.isGeneratingModule(s);
	}

	@Override
	public List<Vector<ComplexNumber>> getModuleGenerators() {
		return vectorSpace.getModuleGenerators();
	}

	@Override
	public Vector<ComplexNumber> asVector(Vector<ComplexNumber> s) {
		return s;
	}

//	@Override
//	public List<Vector<ComplexNumber>> latticeReduction(List<Vector<ComplexNumber>> s) {
//		Rationals q = Rationals.q();
//		NumberField nf = new NumberField(q.getUnivariatePolynomialRing().getPolynomial(q.one(), q.zero(), q.one()));
//		return latticeReduction(s, nf.complexEmbeddings().get(0));
//	}
//
//	public List<Vector<ComplexNumber>> latticeReduction(List<Vector<ComplexNumber>> s, EmbeddedNumberField<ComplexNumber> quadraticField) {
//		return latticeReduction(this, new MathMap<>() {
//			@Override
//			public ComplexNumber evaluate(ComplexNumber t) {
//				return quadraticField.embedding(quadraticField.roundToInteger(t));
//			}
//		}, s);
//	}

	@Override
	public MatrixAlgebra<ComplexNumber> matrixAlgebra() {
		return vectorSpace.matrixAlgebra();
	}

	@Override
	public int dimension() {
		return dimension;
	}

	@Override
	public List<List<ComplexNumber>> nonTrivialCombinations(List<Vector<ComplexNumber>> s) {
		return vectorSpace.nonTrivialCombinations(s);
	}

}
