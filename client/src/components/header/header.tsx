import React, { useContext, useEffect } from "react";
import styles from "./header.module.scss";
import { Link } from "react-router-dom";
import Button from "../../atoms/button/button";
import { AuthService } from "../../services/auth/auth.service";
import Icon from "../../atoms/icons/icon";
import { AppContext } from "../../app-context";
import Logo from "./../../assets/logo.png";

type HeaderProps = { authService: AuthService };

const Header: React.FunctionComponent<HeaderProps> = (props: HeaderProps) => {
  const { isUserLoggedIn, setIsUserLoggedIn } = useContext(AppContext);

  console.log(isUserLoggedIn);
  return (
    <div className={styles["header-container"]}>
      <div className={styles["header-content"]}>
        <div className={styles.btns}>
          <Link className={styles.logo} to="/">
            <img src={Logo} className={styles.logo} />
          </Link>

          {/** VUS */}
          {isUserLoggedIn && (
            <>
              <div className={styles["option-container"]}>
                <div className={styles["option-btn"]}>VUS</div>
                <div className={styles.options}>
                  <Link to={"/vus-upload"} className={styles.option}>
                    VUS Upload
                  </Link>
                  <Link to={"/file-upload"} className={styles.option}>
                    File Upload
                  </Link>
                  <Link to={"/view-vus"} className={styles.option}>
                    View All
                  </Link>
                </div>
              </div>

              {/** Samples */}
              <div className={styles["option-container"]}>
                <div className={styles["option-btn"]}>Samples</div>
                <div className={styles.options}>
                  <Link to={"/view-samples"} className={styles.option}>
                    View All
                  </Link>
                </div>
              </div>
            </>
          )}
        </div>

        {isUserLoggedIn && (
          <div className={styles["side-btns"]}>
            <Link
              to="/login"
              className={styles["logout-btn"]}
              onClick={() =>
                props.authService.logout().then(() => setIsUserLoggedIn(false))
              }
            >
              <Button text="Log out" icon="logout" />
            </Link>
            <Link className={styles.icon} to="/profile">
              <Icon name="profile" fill="#fff" />
            </Link>
          </div>
        )}
      </div>
    </div>
  );
};

export default Header;
