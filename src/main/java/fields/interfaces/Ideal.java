package fields.interfaces;

import java.util.List;

import fields.local.Value;

public interface Ideal<T extends Element<T>> extends Module<T, T>, Element<Ideal<T>> {
	public boolean isPrime();
	public boolean isMaximal();
	public boolean isRadical();
	public List<T> generators();
	public List<T> generate(T t);
	public T residue(T t);
	public boolean contains(T t);
	public Value maximumPowerContains(T t);
	public boolean contains(Ideal<T> t);
	public boolean equalsIdeal(Ideal<T> other);
	public boolean isLeftIdeal();
	public boolean isRightIdeal();
}
