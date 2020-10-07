package fields.interfaces;

import java.util.List;

public interface VectorSpace<T extends Element<T>, S extends Element<S>> extends Module<T, S> {
public Field<T> getField();
public List<S> getBasis();
}
