from flask import current_app
from flask_login import login_user
from werkzeug.security import generate_password_hash, check_password_hash

from server import db
from server.models import ScientificMembers
from server.responses.internal_response import InternalResponse


def signup_scientific_member(email: str, name: str, surname: str, password: str) -> InternalResponse:
    current_app.logger.info(f'Sign up in progress...')

    # if this returns a scientific member, then the email already exists in database
    scientific_member = db.session.query(ScientificMembers).filter(ScientificMembers.email == email).first()

    # if a scientific member is found, we want to redirect back to signup page so user can try again
    if scientific_member:
        current_app.logger.info(f'Scientific member already exists in database!')
        return InternalResponse({'scientificMemberAlreadyExists': True}, 200)

    # create a new scientific member & hash the password
    new_scientific_member = ScientificMembers(email=email, name=name, surname=surname, password=generate_password_hash(password, method='pbkdf2:sha256'))

    # add the new scientific member to the database
    db.session.add(new_scientific_member)
    db.session.commit()

    current_app.logger.info(f'New Scientific member added to database!')
    return InternalResponse({'scientificMemberAlreadyExists': False}, 200)


def login_scientific_member(email: str, password: str, remember: bool) -> InternalResponse:
    current_app.logger.info(f'Login in progress...')

    scientific_member = db.session.query(ScientificMembers).filter(ScientificMembers.email == email).first()

    # check if the user actually exists
    # take the user-supplied password, hash it, and compare it to the hashed password in the database
    if not scientific_member or not check_password_hash(scientific_member.password, password):
        # flash('Please check your login details and try again.')
        current_app.logger.info(f'Login details are not valid!')
        return InternalResponse({'areCredentialsCorrect': False}, 500)

    # if the above check passes, then we know the user has the right credentials
    current_app.logger.info(f'Login details are valid! Logging in user')
    login_user(scientific_member, remember=remember)
    return InternalResponse({'areCredentialsCorrect': True}, 200)