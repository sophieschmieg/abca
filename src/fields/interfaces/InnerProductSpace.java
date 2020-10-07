package fields.interfaces;

import java.util.List;

public interface InnerProductSpace<T extends Element<T>, S extends Element<S>> extends NormedVectorSpace<T, S> {
	T innerProduct(S s1, S s2);
	List<S> gramSchmidt(List<S> s);
	List<S> normedGramSchmidt(List<S> s);
}
