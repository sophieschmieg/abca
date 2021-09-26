package varieties.projective;

import fields.interfaces.Element;
import varieties.Scheme;

public interface ProjectiveScheme<T extends Element<T>> extends Scheme<T, ProjectivePoint<T>> {
	GenericProjectiveScheme<T> asGenericProjectiveScheme();
}
