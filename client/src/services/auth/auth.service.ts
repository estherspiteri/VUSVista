import {
  ILoginRequest,
  ILoginResponse,
  ILogoutResponse,
  IRegisterRequest,
  IRegisterResponse,
} from "./auth.dto";

export class AuthService {
  async isUserLoggedIn(): Promise<ILoginResponse> {
    const result: ILoginResponse = await fetch(`/auth/logged-in-check`, {
      method: "GET",
      headers: {
        "Content-Type": "application/json;charset=UTF-8",
      },
      cache: "no-store",
    })
      .then((response: Response) => {
        return response.json();
      })
      .catch((error) => console.error("error============:", error)); //TODO: handle error

    return result;
  }

  async login(input: ILoginRequest): Promise<ILoginResponse> {
    let data = new FormData();

    // Append the JSON string as a blob to the FormData
    data.append("email", input.email);
    data.append("password", input.password);
    data.append("remember", input.remember.toString());

    const result: ILoginResponse = await fetch(`/auth/login`, {
      method: "POST",
      body: data,
      cache: "no-store",
    })
      .then((response: Response) => {
        return response.json();
      })
      .catch((error) => console.error("error============:", error)); //TODO: handle error

    return result;
  }

  async register(input: IRegisterRequest): Promise<IRegisterResponse> {
    let data = new FormData();

    // Append the JSON string as a blob to the FormData
    data.append("name", input.name);
    data.append("surname", input.surname);
    data.append("email", input.email);
    data.append("password", input.password);

    const result: IRegisterResponse = await fetch(`/auth/register`, {
      method: "POST",
      body: data,
      cache: "no-store",
    })
      .then((response: Response) => {
        return response.json();
      })
      .catch((error) => console.error("error============:", error)); //TODO: handle error

    return result;
  }

  async logout(): Promise<ILogoutResponse> {
    const result: ILogoutResponse = await fetch(`/auth/logout`, {
      method: "POST",
      headers: {
        "Content-Type": "application/json;charset=UTF-8",
      },
    })
      .then((response: Response) => {
        return response.json();
      })
      .catch((error) => console.error("error============:", error)); //TODO: handle error

    return result;
  }
}

export const authService = new AuthService();
