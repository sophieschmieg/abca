package fields.helper;

import fields.interfaces.Element;

public abstract class AbstractElement<T> implements Element<T> {
	@SuppressWarnings("unchecked")
	@Override
	public boolean equals(Object o) {
		return this.compareTo((T)o) == 0;
	}
}
