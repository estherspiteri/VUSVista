export interface ILoginRequest {
  email: string;
  password: string;
  remember: boolean;
}

export interface ILoginResponse {
  isUserLoggedIn: boolean;
}

export interface ILogoutResponse {
  isUserLoggedOut: boolean;
}

export interface IRegisterRequest {
  name: string;
  surname: string;
  email: string;
  password: string;
}

export interface IRegisterResponse {
  scientificMemberAlreadyExists: boolean;
}
