package fields;

public interface MathSet<T extends Element> {
	public T getRandomElement();
	public boolean isFinite();
	public int getNumberOfElements() throws InfinityException;
	public Iterable<T> getElements() throws InfinityException;
}
