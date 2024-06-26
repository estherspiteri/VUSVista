import customFetch from "../api/api.service";
import {
  ILoginRequest,
  ILoginResponse,
  ILogoutResponse,
  IRegisterRequest,
  IRegisterResponse,
} from "./auth.dto";

export class AuthService {
  async isUserLoggedIn(): Promise<ILoginResponse> {
    const result: ILoginResponse = await customFetch(`/auth/logged-in-check`, {
      method: "GET",
      headers: {
        "Content-Type": "application/json;charset=UTF-8",
      },
      cache: "no-store",
    })
      .then((response) => {
        return response
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

    const result: ILoginResponse = await customFetch(`/auth/login`, {
      method: "POST",
      body: data,
      cache: "no-store",
    })
      .then((response) => {
        return response
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

    const result: IRegisterResponse = await customFetch(`/auth/register`, {
      method: "POST",
      body: data,
      cache: "no-store",
    })
      .then((response) => {
        return response
      })
      .catch((error) => console.error("error============:", error)); //TODO: handle error

    return result;
  }

  async logout(): Promise<ILogoutResponse> {
    const result: ILogoutResponse = await customFetch(`/auth/logout`, {
      method: "POST",
      headers: {
        "Content-Type": "application/json;charset=UTF-8",
      },
    })
      .then((response) => {
        return response
      })
      .catch((error) => console.error("error============:", error)); //TODO: handle error

    return result;
  }
}

export const authService = new AuthService();
