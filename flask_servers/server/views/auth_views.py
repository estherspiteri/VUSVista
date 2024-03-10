import json

from flask import Blueprint, request, Response, current_app
from flask_login import login_required, logout_user, current_user

from server.services.auth_service import register_scientific_member, login_scientific_member

auth_views = Blueprint('auth_views', __name__)


@auth_views.route('/logged-in-check', methods=['GET'])
def is_user_logged_in():
    print(current_user.is_authenticated)
    return Response(json.dumps({'isUserLoggedIn': current_user.is_authenticated}), 200, mimetype='application/json')


@auth_views.route('/login', methods=['POST'])
def login():
    email = request.form['email']
    password = request.form['password']
    remember = True if request.form['remember'] == "true" else False

    login_res = login_scientific_member(email, password, remember)
    print(email, password)

    return Response(json.dumps({'isUserLoggedIn': login_res.data['areCredentialsCorrect']}), 200, mimetype='application/json')


@auth_views.route('/register', methods=['POST'])
def register():
    email = request.form['email']
    name = request.form['name']
    surname = request.form['surname']
    password = request.form['password']

    register_res = register_scientific_member(email, name, surname, password)

    return Response(json.dumps({'scientificMemberAlreadyExists': register_res.data['scientificMemberAlreadyExists']}), 200, mimetype='application/json')


@auth_views.route('/logout', methods=['POST'])
@login_required
def logout():
    current_app.logger.info(f"Logging out current user..")

    logout_user()
    return Response(json.dumps({'isUserLoggedOut': True}), 200, mimetype='application/json')

