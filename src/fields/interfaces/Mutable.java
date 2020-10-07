package fields.interfaces;

public interface Mutable<T extends Element<T>> {
	T get();
	boolean isUnmutable();
	T makeUnmutable();
	T cloneUnmutable();
	Mutable<T> clone();
}
