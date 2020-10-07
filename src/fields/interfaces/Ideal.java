package fields.interfaces;

import java.util.List;

public interface Ideal<T extends Element<T>> extends Module<T, T> {
	public boolean isPrime();
	public boolean isMaximal();
	public List<T> generators();
	public List<T> generate(T t);
	public T residue(T t);
	public boolean contains(T t);
	public boolean contains(Ideal<T> t);
	public boolean equalsIdeal(Ideal<T> other);
}
