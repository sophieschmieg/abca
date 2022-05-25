package varieties.projective;

import fields.interfaces.Element;
import fields.interfaces.Field;
import varieties.FunctionField;
import varieties.Scheme;

public interface ProjectiveScheme<T extends Element<T>> extends Scheme<T, ProjectivePoint<T>> {
	Field<T> getField();
	FunctionField<T> getFunctionField();
	GenericProjectiveScheme<T> asGenericProjectiveScheme();
}
