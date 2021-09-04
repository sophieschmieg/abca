package fields.interfaces;

import java.util.List;

public interface Algebra<T extends Element<T>, S extends Element<S>> extends Module<T, S>, Ring<S> {
	public S getEmbedding(T t);
	
	public boolean isGeneratingAlgebra(List<S> s);

	public List<S> getAlgebraGenerators();

}
