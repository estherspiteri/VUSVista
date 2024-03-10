import React, { useState } from "react";
import styles from "./login-page.module.scss";
import Button from "../../atoms/button/button";
import { AuthService } from "../../services/auth/auth.service";
import Text from "../../atoms/text/text";
import { Link } from "react-router-dom";

type LoginPageProps = { authService: AuthService };

const LoginPage: React.FunctionComponent<LoginPageProps> = (
  props: LoginPageProps
) => {
  const [email, setEmail] = useState("");
  const [emailErrorMsg, setEmailErrorMsg] = useState("");

  const [password, setPassword] = useState("");
  const [passwordErrorMsg, setPasswordErrorMsg] = useState("");

  const [loginErrorMsg, setLoginErrorMsg] = useState("");

  const [rememberMe, setRememberMe] = useState(false);

  return (
    <div className={styles["login-page-container"]}>
      <div className={styles.title}>Login</div>

      <div className={styles["login-content"]}>
        <div className={styles["field-container"]}>
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

          <span className={styles["remember-me"]}>
            <input
              type="checkbox"
              defaultChecked={rememberMe}
              name="rememberMe"
              onChange={() => setRememberMe(!rememberMe)}
            />
            Remember me
          </span>
        </div>

        <div className={styles["login-btn-container"]}>
          {loginErrorMsg.length > 0 && (
            <p className={styles.error}>{loginErrorMsg}</p>
          )}
          <Button text="Login" onClick={login} className={styles.btn} />
        </div>
      </div>

      <p className={styles.register}>
        <Link to={"/register"} className={styles.link}>
          Click here
        </Link>
        &nbsp;to register.
      </p>
    </div>
  );

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
      setPasswordErrorMsg("Please enter your email");
    }

    return isValid;
  }

  function validateCredentials() {
    const isEmailValid = checkEmailValidity(email);
    const isPasswordValid = checkPasswordValidity(password);

    return isEmailValid && isPasswordValid;
  }

  function login() {
    const isValid = validateCredentials();

    if (isValid) {
      props.authService
        .login({ email: email, password: password, remember: rememberMe })
        .then((res) => {
          console.log(res);
          if (res.isUserLoggedIn) {
            //redirect to profile page
          } else {
            setLoginErrorMsg(
              "The credentials you have inputted do not match any existing user. Kindly check and try again."
            );
          }
        });
    }
  }
};

export default LoginPage;
