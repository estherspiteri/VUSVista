import React, { useState } from "react";
import styles from "./header.module.scss";
import { Link } from "react-router-dom";
import Button from "../../atoms/button/button";
import { AuthService } from "../../services/auth/auth.service";

type HeaderProps = { isUserLoggedIn: boolean; authService: AuthService };

const Header: React.FunctionComponent<HeaderProps> = (props: HeaderProps) => {
  return (
    <div className={styles["header-container"]}>
      <div className={styles["header-content"]}>
        {props.isUserLoggedIn && (
          <>
            <div className={styles.btns}>
              {/** VUS */}
              <div className={styles["option-container"]}>
                <div className={styles["option-btn"]}>VUS</div>
                <div className={styles.options}>
                  <Link to={"/view-vus"} className={styles.option}>
                    View All
                  </Link>
                </div>
              </div>

              {/** Samples */}
              <div className={styles["option-container"]}>
                <div className={styles["option-btn"]}>Samples</div>
                <div className={styles.options}>
                  <Link to={"upload-vus"} className={styles.option}>
                    File Upload
                  </Link>
                  <Link to={"/view-samples"} className={styles.option}>
                    View All
                  </Link>
                </div>
              </div>
            </div>

            <div className={styles["logout-btn"]}>
              <Button
                text="Log out"
                icon="logout"
                onClick={() => props.authService.logout()} //TODO: redirect to login
              />
            </div>
          </>
        )}
      </div>
    </div>
  );
};

export default Header;
