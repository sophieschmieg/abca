package fields.vectors;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;

import fields.floatingpoint.Complex.ComplexNumber;
import fields.floatingpoint.Reals.Real;
import fields.helper.AbstractElement;
import fields.interfaces.Element;
import fields.interfaces.InnerProductSpace;
import fields.interfaces.MathMap;
import fields.interfaces.ValueField;
import fields.vectors.DualVectorSpace.Dual;

public class DualVectorSpace<T extends Element<T>, S extends Element<S>>
		extends AbstractInnerProductSpace<T, Dual<T, S>> implements InnerProductSpace<T, Dual<T, S>> {
	public static class Dual<T extends Element<T>, S extends Element<S>> extends AbstractElement<Dual<T, S>>
			implements MathMap<S, T> {
		private InnerProductSpace<T, S> space;
		private S dual;

		private Dual(InnerProductSpace<T, S> space, S dual) {
			this.space = space;
			this.dual = dual;
		}

		@Override
		public int compareTo(Dual<T, S> o) {
			return dual.compareTo(o.dual);
		}

		@Override
		public T evaluate(S t) {
			return space.innerProduct(dual, t);
		}

		public S dual() {
			return dual;
		}

		public Matrix<T> asMatrix() {
			return Matrix.fromRows(Collections.singletonList(space.asVector(dual)));
		}

		public String toString() {
			return dual.toString() + "^t";
		}
	}

	private InnerProductSpace<T, S> space;

	public DualVectorSpace(InnerProductSpace<T, S> space) {
		this.space = space;
	}

	public String toString() {
		return space.toString() + "'";
	}

	@Override
	public ComplexNumber asComplexNumber(T t) {
		return space.asComplexNumber(t);
	}

	@Override
	public T fromComplexNumber(ComplexNumber t) {
		return space.fromComplexNumber(t);
	}

	@Override
	public T fromReal(Real r) {
		return space.fromReal(r);
	}

	@Override
	public ValueField<T> getValueField() {
		return space.getValueField();
	}

	@Override
	public DualVectorSpace<T, S> withDimension(int dimension) {
		return new DualVectorSpace<>(space.withDimension(dimension));
	}

	public Dual<T, S> canonicalIsomorphism(S s) {
		return new Dual<>(space, s);
	}

	public Dual<T, S> fromMatrix(Matrix<T> s) {
		if (s.rows() != 1) {
			throw new ArithmeticException("Expected only one row!");
		}
		return canonicalIsomorphism(/* space.conjugate( */space.fromVector(s.row(1)));
	}

	public Dual<T, S> fromRowVector(Vector<T> s) {
		return fromMatrix(Matrix.fromRows(Collections.singletonList(s)));
	}

	@Override
	public T conjugate(T s) {
		return space.conjugate(s);
	}

	@Override
	public Dual<T, S> getUnitVector(int index) {
		return canonicalIsomorphism(space.getUnitVector(index));
	}

	public List<Dual<T, S>> getDualBasis(List<S> basis) {
		List<Vector<T>> asColumns = new ArrayList<>();
		for (S basisVector : basis) {
			asColumns.add(space.asVector(basisVector));
		}
		Matrix<T> baseChangeMatrix = Matrix.fromColumns(asColumns);
		return getDualBasis(baseChangeMatrix);
	}

	public List<Dual<T, S>> getDualBasis(Matrix<T> baseChangeMatrix) {
		Matrix<T> inverseBaseChange = space.matrixAlgebra().inverse(baseChangeMatrix);
		List<Dual<T, S>> result = new ArrayList<>();
		for (int i = 0; i < dimension(); i++) {
			result.add(fromRowVector(inverseBaseChange.row(i + 1)));
		}
		return result;
	}

	@Override
	public List<Dual<T, S>> getBasis() {
		List<Dual<T, S>> result = new ArrayList<>();
		for (S b : space.getBasis()) {
			result.add(canonicalIsomorphism(b));
		}
		return result;
	}

	@Override
	public int dimension() {
		return space.dimension();
	}

	@Override
	public Vector<T> asVector(Dual<T, S> s) {
		return space.asVector(s.dual);
	}

	@Override
	public MatrixAlgebra<T> matrixAlgebra() {
		return space.matrixAlgebra();
	}

	@Override
	public Dual<T, S> zero() {
		return canonicalIsomorphism(space.zero());
	}

	@Override
	public Dual<T, S> add(Dual<T, S> s1, Dual<T, S> s2) {
		return canonicalIsomorphism(space.add(s1.dual, s2.dual));
	}

	@Override
	public Dual<T, S> negative(Dual<T, S> s) {
		return canonicalIsomorphism(space.negative(s.dual));
	}

	@Override
	public Dual<T, S> scalarMultiply(T t, Dual<T, S> s) {
		return canonicalIsomorphism(space.scalarMultiply(t, s.dual));
	}

	private List<S> dualList(List<Dual<T, S>> list) {
		List<S> result = new ArrayList<>();
		for (Dual<T, S> v : list) {
			result.add(v.dual);
		}
		return result;
	}

	@Override
	public boolean isLinearIndependent(List<Dual<T, S>> s) {
		return space.isLinearIndependent(dualList(s));
	}

	@Override
	public boolean isGeneratingModule(List<Dual<T, S>> s) {
		return space.isGeneratingModule(dualList(s));
	}

	@Override
	public List<List<T>> nonTrivialCombinations(List<Dual<T, S>> s) {
		return space.nonTrivialCombinations(dualList(s));
	}

	@Override
	public Exactness exactness() {
		return space.exactness();
	}

	@Override
	public Dual<T, S> getRandomElement() {
		return canonicalIsomorphism(space.getRandomElement());
	}

	@Override
	public Iterator<Dual<T, S>> iterator() {
		return new Iterator<>() {
			private Iterator<S> it = space.iterator();

			@Override
			public boolean hasNext() {
				return it.hasNext();
			}

			@Override
			public Dual<T, S> next() {
				return canonicalIsomorphism(it.next());
			}
		};
	}

	@Override
	public T innerProduct(Dual<T, S> s1, Dual<T, S> s2) {
		return space.innerProduct(s2.dual, s1.dual);
	}
}
