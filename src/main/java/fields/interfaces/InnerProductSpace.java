package fields.interfaces;

import java.util.ArrayList;
import java.util.List;

import fields.floatingpoint.Reals;
import fields.floatingpoint.Reals.Real;

public interface InnerProductSpace<T extends Element<T>, S extends Element<S>> extends NormedVectorSpace<T, S> {
	T innerProduct(S s1, S s2);

	@Override
	default Real valueNorm(S s) {
		Reals r = Reals.r(1024);
		return r.positiveSqrt(getValueField().value(innerProduct(s, s)));
	}

	default List<S> gramSchmidt(List<S> s) {
		List<S> orthogonal = new ArrayList<>();
		ValueField<T> f = getValueField();
		for (int i = 0; i < s.size(); i++) {
			S vector = s.get(i);
			for (int j = 0; j < i; j++) {
				vector = subtract(vector, scalarMultiply(
						f.divide(innerProduct(s.get(i), s.get(j)), innerProduct(s.get(j), s.get(j))), s.get(j)));
			}
			orthogonal.add(vector);
		}
		return orthogonal;

	}

	default List<S> normedGramSchmidt(List<S> s) {
		List<S> orthogonal = gramSchmidt(s);
		List<S> orthonormal = new ArrayList<>();
		ValueField<T> f = getValueField();
		for (S vector : orthogonal) {
			orthonormal.add(
					scalarMultiply(f.inverse(f.sqrt(innerProduct(vector, vector)).keySet().iterator().next()), vector));
		}
		return orthonormal;
	}

	List<S> latticeReduction(List<S> s);
}
