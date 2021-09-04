package fields.exceptions;

public class UnmodifiableElementException extends RuntimeException {
	private static final long serialVersionUID = 1L;

	public UnmodifiableElementException() {
		super("Element is not modifiable");
	}
}
