import React, { useState } from "react";
import styles from "./register-page.module.scss";
import Button from "../../atoms/button/button";
import { AuthService } from "../../services/auth/auth.service";
import Text from "../../atoms/text/text";
import { Link } from "react-router-dom";

type RegisterPageProps = { authService: AuthService };

const RegisterPage: React.FunctionComponent<RegisterPageProps> = (
  props: RegisterPageProps
) => {
  const [name, setName] = useState("");
  const [nameErrorMsg, setNameErrorMsg] = useState("");

  const [surname, setSurname] = useState("");
  const [surnameErrorMsg, setSurnameErrorMsg] = useState("");

  const [email, setEmail] = useState("");
  const [emailErrorMsg, setEmailErrorMsg] = useState("");

  const [password, setPassword] = useState("");
  const [passwordErrorMsg, setPasswordErrorMsg] = useState("");

  const [registerErrorMsg, setRegisterErrorMsg] = useState("");

  return (
    <div className={styles["register-page-container"]}>
      <div className={styles.title}>Register</div>

      <div className={styles["register-content"]}>
        <div className={styles["field-container"]}>
          {/** Name */}
          <div className={styles.field}>
            <span className={styles["field-name"]}>Name:</span>
            <Text
              value={name}
              errorMsg={nameErrorMsg}
              onChange={(e) => setName(e.currentTarget.value)}
              validationCallback={(val) => checkNameValidity(val, true)}
            />
          </div>

          {/** Surname */}
          <div className={styles.field}>
            <span className={styles["field-name"]}>Surname:</span>
            <Text
              value={surname}
              errorMsg={surnameErrorMsg}
              onChange={(e) => setSurname(e.currentTarget.value)}
              validationCallback={(val) => checkNameValidity(val, false)}
            />
          </div>

          {/** Email - TODO: add validation */}
          <div className={styles.field}>
            <span className={styles["field-name"]}>Email:</span>
            <Text
              value={email}
              errorMsg={emailErrorMsg}
              onChange={(e) => setEmail(e.currentTarget.value)}
              validationCallback={checkEmailValidity}
            />
          </div>

          {/** Password - TODO: add visibility */}
          <div className={styles.field}>
            <span className={styles["field-name"]}>Password:</span>
            <Text
              type="password"
              value={password}
              errorMsg={passwordErrorMsg}
              onChange={(e) => setPassword(e.currentTarget.value)}
              validationCallback={checkPasswordValidity}
            />
          </div>
        </div>

        <div className={styles["register-btn-container"]}>
          {registerErrorMsg.length > 0 && (
            <p className={styles.error}>{registerErrorMsg}</p>
          )}
          <Button text="Register" onClick={register} className={styles.btn} />
        </div>
      </div>

      <p className={styles.login}>
        <Link to={"/login"} className={styles.link}>
          Click here
        </Link>
        &nbsp;to login.
      </p>
    </div>
  );

  function checkNameValidity(val: string, isFirstName: boolean) {
    const isValid = val.trim().length > 0;
    let errorMsg = "";

    if (!isValid) {
      errorMsg = `Please enter your ${isFirstName ? "name" : "surname"}`;
    }

    if (isFirstName) {
      setNameErrorMsg(errorMsg);
    } else {
      setSurnameErrorMsg(errorMsg);
    }

    return isValid;
  }

  function checkEmailValidity(val: string) {
    const isValid = val.trim().length > 0;

    if (isValid) {
      setEmailErrorMsg("");
    } else {
      setEmailErrorMsg("Please enter your email");
    }

    return isValid;
  }

  function checkPasswordValidity(val: string) {
    const isValid = val.trim().length > 0;

    if (isValid) {
      setPasswordErrorMsg("");
    } else {
      setPasswordErrorMsg("Please enter your password");
    }

    return isValid;
  }

  function validateCredentials() {
    const isNameValid = checkNameValidity(name, true);
    const isSurnameValid = checkNameValidity(surname, false);
    const isEmailValid = checkEmailValidity(email);
    const isPasswordValid = checkPasswordValidity(password);

    return isNameValid && isSurnameValid && isEmailValid && isPasswordValid;
  }

  function register() {
    const isValid = validateCredentials();

    if (isValid) {
      props.authService
        .register({
          name: name,
          surname: surname,
          email: email,
          password: password,
        })
        .then((res) => {
          if (res.scientificMemberAlreadyExists) {
            setRegisterErrorMsg(
              "A user with the inputted email already exists."
            );
          } else {
            window.location.href = "/login";
          }
        });
    }
  }
};

export default RegisterPage;
